using Printf
using LinearAlgebra
using JuMP
import MathOptInterface as MOI
import StaticArrays: SVector
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim
const SC = OP.SymbolicCertifier
const SY = DI.Symbolic

module PBRef
include(joinpath(@__DIR__, "..", "Permis_B", "permis_B.jl"))
end

module PERef
include(joinpath(@__DIR__, "get_nominal_traje.jl"))
end

const UNWRAP_PE = true
const LMI_CFG_REF = Ref{Any}(nothing)
const OPTS_REF = Ref{Any}(nothing)
const WDOM_REF = Ref{Any}(nothing)

function make_problem_from_pb_run(pb_run)::PR.OptimalControlProblem
    ctrl = PBRef.build_control_problem()
    return PR.OptimalControlProblem(
        pb_run.concrete_system,
        ctrl.initial_set,
        ctrl.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

function candidate_from_lists(xs, us, Ts; source)
    return OP.CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = Ts,
        source = source,
        metadata = (;),
    )
end

function _make_lmi_cfg_single_thread(base_cfg)
    backend = base_cfg.sdp_opt
    if occursin("Mosek", string(typeof(backend))) && isdefined(PBRef, :MosekTools)
        silent = !base_cfg.verbose
        log_level = base_cfg.verbose ? 10 : 0
        backend = optimizer_with_attributes(
            PBRef.MosekTools.Optimizer,
            MOI.Silent() => silent,
            MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => log_level,
            MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
        )
        println("[info] MOSEK backend forced to single thread (MSK_IPAR_NUM_THREADS=1)")
    else
        @warn "Could not enforce MOSEK thread option on lmi_cfg.sdp_opt; backend may be non-MOSEK or inaccessible."
    end

    return PBRef.LMIConfig(
        rayon_terminal = base_cfg.rayon_terminal,
        λ = base_cfg.λ,
        maxδx = base_cfg.maxδx,
        maxδu = base_cfg.maxδu,
        sdp_opt = backend,
        symbolic_rk4_substeps = base_cfg.symbolic_rk4_substeps,
        ΔX = base_cfg.ΔX,
        ΔU = base_cfg.ΔU,
        ΔW = base_cfg.ΔW,
        verbose = base_cfg.verbose,
    )
end

function run_new(xs, us, Ts, problem, av_mod, params)
    lmi_cfg = LMI_CFG_REF[]
    opts = OPTS_REF[]
    Wdom = WDOM_REF[]

    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        lmi_cfg.sdp_opt,
        opts,
    )

    cand = candidate_from_lists(xs, us, Ts; source = :quad_compare)

    builder = function (prob, candidate, cfg)
        o = cfg.options
        return av_mod.symbolic_system(
            prob.system.X;
            _U_ = prob.system.U,
            params = params,
            Ts = candidate.Ts,
            ΔX = o.ΔX,
            ΔU = o.ΔU,
            ΔW = o.ΔW,
            rk4_num_substeps = o.symbolic_rk4_substeps,
        )
    end

    ctx = SC.build_symbolic_context(problem, cand, cert_cfg, builder)
    return SC.run_backward_chain!(ctx)
end

function run_old(xs, us, Ts, problem, av_mod, params)
    lmi_cfg = LMI_CFG_REF[]

    run_like = (;
        concrete_system = problem.system,
        params = params,
        Δt = Ts,
        state_list = xs,
        input_list = us,
    )

    old_ctx = PBRef.build_symbolic_lmi_context(run_like, lmi_cfg)

    E_terminal = if isdefined(PBRef, :creer_ellipsoide_terminal)
        PBRef.creer_ellipsoide_terminal(xs; rayon = lmi_cfg.rayon_terminal)
    else
        nx = length(xs[end])
        P = Matrix{Float64}(I, nx, nx) * (1.0 / lmi_cfg.rayon_terminal^2)
        UT.Ellipsoid(P, collect(Float64, xs[end]))
    end

    return PBRef.synthetiser_transitions_backward(
        xs,
        us,
        old_ctx;
        E_terminal = E_terminal,
        verbose = false,
    )
end

function summarize_old(transitions)
    ells = transitions.ellipsoides
    vol_terminal = isempty(ells) ? NaN : UT.get_volume(ells[1])
    vol_initial = isempty(ells) ? NaN : UT.get_volume(ells[end])
    return (;
        success = Bool(transitions.success),
        failed_k = transitions.failed_k,
        n_ok_steps = length(transitions.indices),
        vol_terminal,
        vol_initial,
    )
end

function summarize_new(res)
    steps = res.steps
    n_ok_steps = count(s -> s.status == :ok, steps)

    vol_terminal = NaN
    vol_initial = NaN
    if res.lmi_data !== nothing && hasproperty(res.lmi_data, :ellipsoids)
        ells = getproperty(res.lmi_data, :ellipsoids)
        if !isempty(ells)
            vol_terminal = UT.get_volume(ells[1])
            vol_initial = UT.get_volume(ells[end])
        end
    end

    return (;
        success = Bool(res.success),
        failed_k = res.failed_k,
        n_ok_steps,
        vol_terminal,
        vol_initial,
    )
end

function _fmt_failed_k(x)
    return x === nothing ? "-" : string(x)
end

function _fmt_float(x)
    return isnan(x) ? "NaN" : @sprintf("%.6e", x)
end

function _run_case(label, traj_name, cert_name, f)
    t0 = time()
    try
        out = f()
        elapsed = time() - t0
        return (;
            label,
            traj = traj_name,
            cert = cert_name,
            ok = true,
            out,
            elapsed,
            err = nothing,
        )
    catch err
        elapsed = time() - t0
        return (;
            label,
            traj = traj_name,
            cert = cert_name,
            ok = false,
            out = nothing,
            elapsed,
            err,
        )
    end
end

function main()
    BLAS.set_num_threads(1)

    # 1) Shared LMI config
    lmi_cfg = _make_lmi_cfg_single_thread(PBRef.LMIConfig())
    opts = (
        maxδx = lmi_cfg.maxδx,
        maxδu = lmi_cfg.maxδu,
        λ = lmi_cfg.λ,
        rayon_terminal = lmi_cfg.rayon_terminal,
        ΔX = lmi_cfg.ΔX,
        ΔU = lmi_cfg.ΔU,
        ΔW = lmi_cfg.ΔW,
        symbolic_rk4_substeps = lmi_cfg.symbolic_rk4_substeps,
    )
    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))

    LMI_CFG_REF[] = lmi_cfg
    OPTS_REF[] = opts
    WDOM_REF[] = Wdom

    println("=== Quad Compare Trajectory x Certifier ===")
    println("LMI cfg: λ=", lmi_cfg.λ, " | maxδx=", lmi_cfg.maxδx, " | maxδu=", lmi_cfg.maxδu)
    println("ΔW=", lmi_cfg.ΔW)

    # 2) Load Permis_B trajectory
    pb_run = PBRef.run_nominal_simulation(;
        use_cache = true,
        force_recompute = false,
        save_plots = false,
        show_animation = false,
    )
    pb_lmi = PBRef.prepare_run_result_for_lmi(pb_run; enable_unwrap = true)
    xs_pb = pb_lmi.state_list
    us_pb = pb_lmi.input_list
    Ts_pb = pb_run.Δt
    problem_pb = make_problem_from_pb_run(pb_run)
    params_pb = pb_run.params
    av_mod_pb = PBRef.AV

    # 3) Load Permis_E trajectory
    pe_nom = PERef.get_nominal_traje()
    @assert pe_nom.candidate !== nothing
    xs_pe_raw = collect(ST.enum_elems(pe_nom.candidate.x_traj))
    us_pe = collect(ST.enum_elems(pe_nom.candidate.u_traj))
    xs_pe = UNWRAP_PE ? PBRef.unwrap_periodic_state_list(xs_pe_raw, SVector(3, 4), SVector(2pi, 2pi)) : xs_pe_raw
    Ts_pe = pe_nom.candidate.Ts
    problem_pe = pe_nom.problem

    av_mod_pe = isdefined(PERef, :AV) ? PERef.AV : PBRef.AV
    params_pe_new = if hasproperty(pe_nom, :params)
        getproperty(pe_nom, :params)
    else
        av_mod_pe.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    end
    # old helper path requires PBRef.AV.Params type
    params_pe_old = PBRef.AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)

    println("Permis_B traj: n_states=", length(xs_pb), " | n_inputs=", length(us_pb))
    println("Permis_E traj: n_states=", length(xs_pe), " | n_inputs=", length(us_pe), " | unwrap=", UNWRAP_PE)

    # 4) Run 2x2 matrix
    c1 = _run_case("PB-old", "Permis_B", "old") do
        run_old(xs_pb, us_pb, Ts_pb, problem_pb, av_mod_pb, params_pb)
    end
    c2 = _run_case("PB-new", "Permis_B", "new") do
        run_new(xs_pb, us_pb, Ts_pb, problem_pb, av_mod_pb, params_pb)
    end
    c3 = _run_case("PE-old", "Permis_E", "old") do
        run_old(xs_pe, us_pe, Ts_pe, problem_pe, av_mod_pe, params_pe_old)
    end
    c4 = _run_case("PE-new", "Permis_E", "new") do
        run_new(xs_pe, us_pe, Ts_pe, problem_pe, av_mod_pe, params_pe_new)
    end

    rows = [c1, c2, c3, c4]

    println("\n=== Results ===")
    @printf("%-10s | %-8s | %-4s | %-7s | %-8s | %-9s | %-13s | %-13s | %-8s\n",
        "Case",
        "Traj",
        "Cert",
        "ran_ok",
        "success",
        "failed_k",
        "n_ok_steps",
        "vol_terminal",
        "vol_init",
    )
    println(repeat("-", 104))

    for r in rows
        if !r.ok
            @printf("%-10s | %-8s | %-4s | %-7s | %-8s | %-9s | %-13s | %-13s | %-8s\n",
                r.label,
                r.traj,
                r.cert,
                "false",
                "ERR",
                "-",
                "-",
                "-",
                "-",
            )
            println("  error(", r.label, "): ", sprint(showerror, r.err))
            continue
        end

        s = r.cert == "old" ? summarize_old(r.out) : summarize_new(r.out)
        @printf("%-10s | %-8s | %-4s | %-7s | %-8s | %-9s | %-13d | %-13s | %-8s\n",
            r.label,
            r.traj,
            r.cert,
            "true",
            string(s.success),
            _fmt_failed_k(s.failed_k),
            s.n_ok_steps,
            _fmt_float(s.vol_terminal),
            _fmt_float(s.vol_initial),
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(1)
    main()
end
