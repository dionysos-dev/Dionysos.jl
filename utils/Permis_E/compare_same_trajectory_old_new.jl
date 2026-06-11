using Printf
using LinearAlgebra
import IntervalArithmetic as IA
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

Base.@kwdef struct CompareConfig
    traj::Symbol = :permis_e
    nsteps::Int = 0
    unwrap::Bool = true
    verbose::Bool = false
end

function parse_args(args)
    traj = :permis_e
    nsteps = 0
    unwrap = true
    verbose = false
    for a in args
        if startswith(a, "--traj=")
            v = lowercase(split(a, "="; limit = 2)[2])
            if v == "permis_e"
                traj = :permis_e
            elseif v == "permis_b"
                traj = :permis_b
            else
                error("--traj must be permis_e or permis_b")
            end
        elseif startswith(a, "--nsteps=")
            nsteps = parse(Int, split(a, "="; limit = 2)[2])
        elseif a == "--no-unwrap"
            unwrap = false
        elseif a == "--verbose"
            verbose = true
        else
            error("Unknown argument: $a")
        end
    end
    return CompareConfig(; traj, nsteps, unwrap, verbose)
end

_vecf(x) = collect(Float64, x)
_matf(X) = Matrix{Float64}(X)

function _maxabs(v)
    isempty(v) && return 0.0
    return maximum(abs, v)
end

function _maxabsdiff_vec(x, y)
    xv = _vecf(x)
    yv = _vecf(y)
    length(xv) == length(yv) || return Inf
    return _maxabs(xv .- yv)
end

function _maxabsdiff_mat(X, Y)
    Xm = _matf(X)
    Ym = _matf(Y)
    size(Xm) == size(Ym) || return Inf
    return _maxabs(Xm .- Ym)
end

function _make_problem_from_pb_run(run_result)
    ctrl = PBRef.build_control_problem()
    return PR.OptimalControlProblem(
        run_result.concrete_system,
        ctrl.initial_set,
        ctrl.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

function _load_trajectory(cfg::CompareConfig)
    if cfg.traj == :permis_b
        pb_run = PBRef.run_nominal_simulation(;
            use_cache = true,
            force_recompute = false,
            save_plots = false,
            show_animation = false,
        )
        pb_lmi = PBRef.prepare_run_result_for_lmi(pb_run; enable_unwrap = cfg.unwrap)
        xs = pb_lmi.state_list
        us = pb_lmi.input_list
        Ts = pb_run.Δt
        problem = _make_problem_from_pb_run(pb_run)
        return (;
            xs,
            us,
            Ts,
            problem,
            params_pb = pb_run.params,
            params_pe = nothing,
            av_mod = PBRef.AV,
        )
    end

    pe_nom = PERef.get_nominal_traje()
    pe_nom.candidate === nothing && error("Permis_E candidate is nothing.")
    xs_raw = collect(ST.enum_elems(pe_nom.candidate.x_traj))
    us = collect(ST.enum_elems(pe_nom.candidate.u_traj))
    xs =
        cfg.unwrap ?
        PBRef.unwrap_periodic_state_list(xs_raw, SVector(3, 4), SVector(2pi, 2pi)) : xs_raw

    params_pb = PBRef.AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    params_pe = PERef.AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    return (;
        xs,
        us,
        Ts = pe_nom.candidate.Ts,
        problem = pe_nom.problem,
        params_pb,
        params_pe,
        av_mod = PERef.AV,
    )
end

function _make_new_context(problem, xs, us, Ts, av_mod, params, lmi_cfg)
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
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        lmi_cfg.sdp_opt,
        opts,
    )

    cand = OP.CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = Ts,
        source = :compare_same_traj,
        metadata = (;),
    )

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

    return SC.build_symbolic_context(problem, cand, cert_cfg, builder)
end

function _linearize_ctx(ctx, xk, uk)
    wk = zeros(length(ctx.w))
    Xbar = IA.IntervalBox(_vecf(xk) .+ ctx.ΔX)
    Ubar = IA.IntervalBox(_vecf(uk) .+ ctx.ΔU)
    Wbar = IA.IntervalBox(wk .+ ctx.ΔW)

    affineSys, L = ST.buildAffineApproximation(
        ctx.fsymbolic,
        ctx.x,
        ctx.u,
        ctx.w,
        _vecf(xk),
        _vecf(uk),
        wk,
        Xbar,
        Ubar,
        Wbar,
    )
    return affineSys, L
end

_backend(ctx) =
    hasproperty(ctx, :backend) ? getproperty(ctx, :backend) : getproperty(ctx, :sdp_opt)

function _terminal_ellipsoid(xs, radius)
    nx = length(xs[end])
    P = Matrix{Float64}(I, nx, nx) * (1.0 / radius^2)
    return UT.Ellipsoid(P, _vecf(xs[end]))
end

function compare_stepwise(
    old_ctx,
    new_ctx,
    xs,
    us,
    lmi_cfg;
    nsteps::Int = 0,
    verbose::Bool = false,
)
    n_u = length(us)
    n_use = nsteps <= 0 ? n_u : min(nsteps, n_u)

    E_old = _terminal_ellipsoid(xs, lmi_cfg.rayon_terminal)
    E_new = _terminal_ellipsoid(xs, lmi_cfg.rayon_terminal)

    println("\n=== Stepwise comparison on SAME trajectory ===")
    println("n_states=", length(xs), " | n_inputs=", n_u, " | compared_steps=", n_use)

    first_mismatch = nothing
    common_failure = nothing
    for j in 1:n_use
        i_state = n_u - j + 1
        xk = xs[i_state]
        uk = us[i_state]

        aff_old, L_old = _linearize_ctx(old_ctx, xk, uk)
        aff_new, L_new = _linearize_ctx(new_ctx, xk, uk)

        dA = _maxabsdiff_mat(aff_old.A, aff_new.A)
        dB = _maxabsdiff_mat(aff_old.B, aff_new.B)
        dc_aff = _maxabsdiff_vec(aff_old.c, aff_new.c)
        dD = _maxabsdiff_mat(aff_old.D, aff_new.D)
        dL = _maxabsdiff_vec(L_old, L_new)

        E_old_prev, _, cost_old = SY.transition_backward(
            aff_old,
            E_old,
            _vecf(xk),
            _vecf(uk),
            old_ctx.Uformat,
            old_ctx.Wformat,
            old_ctx.S,
            L_old,
            _backend(old_ctx);
            maxδx = old_ctx.maxδx,
            maxδu = old_ctx.maxδu,
            λ = old_ctx.λ,
        )

        E_new_prev, _, cost_new = SY.transition_backward(
            aff_new,
            E_new,
            _vecf(xk),
            _vecf(uk),
            new_ctx.Uformat,
            new_ctx.Wformat,
            new_ctx.S,
            L_new,
            _backend(new_ctx);
            maxδx = new_ctx.maxδx,
            maxδu = new_ctx.maxδu,
            λ = new_ctx.λ,
        )

        st_old = E_old_prev === nothing ? :fail : :ok
        st_new = E_new_prev === nothing ? :fail : :ok

        if st_old == :ok && st_new == :ok
            dP = _maxabsdiff_mat(E_old_prev.P, E_new_prev.P)
            dc_ell = _maxabsdiff_vec(E_old_prev.c, E_new_prev.c)
            dcost = abs(Float64(cost_old) - Float64(cost_new))
            @printf(
                "j=%d i=%d status=ok/ok | dA=%.1e dB=%.1e dc=%.1e dD=%.1e dL=%.1e | dP=%.1e dc_ell=%.1e dcost=%.1e\n",
                j,
                i_state,
                dA,
                dB,
                dc_aff,
                dD,
                dL,
                dP,
                dc_ell,
                dcost,
            )
            E_old = E_old_prev
            E_new = E_new_prev
        else
            @printf(
                "j=%d i=%d status=%s/%s | dA=%.1e dB=%.1e dc=%.1e dD=%.1e dL=%.1e\n",
                j,
                i_state,
                string(st_old),
                string(st_new),
                dA,
                dB,
                dc_aff,
                dD,
                dL,
            )
            if st_old == st_new
                common_failure = (; j, i_state, status = st_old, dA, dB, dc_aff, dD, dL)
            else
                first_mismatch = (; j, i_state, st_old, st_new, dA, dB, dc_aff, dD, dL)
            end
            break
        end

        if verbose
            @printf(
                "  volumes old/new: %.12e / %.12e\n",
                UT.get_volume(E_old),
                UT.get_volume(E_new)
            )
        end
    end

    if first_mismatch !== nothing
        println(
            "Result: OLD/NEW mismatch at j=",
            first_mismatch.j,
            " (i_state=",
            first_mismatch.i_state,
            ")",
        )
    elseif common_failure !== nothing
        println(
            "Result: identical behavior, common failure at j=",
            common_failure.j,
            " (i_state=",
            common_failure.i_state,
            ")",
        )
    else
        println("Result: no mismatch found on tested steps.")
    end
end

function main(args = ARGS)
    cfg = parse_args(args)
    lmi_cfg = PBRef.LMIConfig()

    data = _load_trajectory(cfg)
    xs = data.xs
    us = data.us

    println("=== Compare old/new on identical trajectory ===")
    println("traj source = ", cfg.traj)
    println("unwrap = ", cfg.unwrap)
    println("λ=", lmi_cfg.λ, " | maxδx=", lmi_cfg.maxδx, " | maxδu=", lmi_cfg.maxδu)
    println("ΔW=", lmi_cfg.ΔW)

    run_like = (;
        concrete_system = data.problem.system,
        params = data.params_pb,
        Δt = data.Ts,
        state_list = xs,
        input_list = us,
    )
    old_ctx = PBRef.build_symbolic_lmi_context(run_like, lmi_cfg)

    params_new = data.params_pe === nothing ? data.params_pb : data.params_pe
    new_ctx =
        _make_new_context(data.problem, xs, us, data.Ts, data.av_mod, params_new, lmi_cfg)

    return compare_stepwise(
        old_ctx,
        new_ctx,
        xs,
        us,
        lmi_cfg;
        nsteps = cfg.nsteps,
        verbose = cfg.verbose,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
