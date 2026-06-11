using LinearAlgebra
using Printf
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

function _parse_args(args)
    nsteps = 15
    skip_permis_e = false
    for a in args
        if startswith(a, "--nsteps=")
            nsteps = parse(Int, split(a, "="; limit = 2)[2])
        elseif a == "--skip-permis-e"
            skip_permis_e = true
        else
            error("Unknown argument: $a")
        end
    end
    return (; nsteps, skip_permis_e)
end

_vecf(x) = collect(Float64, x)
_matf(M) = Matrix{Float64}(M)

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

function _trajectory_delta(xs1, us1, xs2, us2)
    nx = min(length(xs1), length(xs2))
    nu = min(length(us1), length(us2))

    dx = 0.0
    ix = 1
    for i in 1:nx
        di = _maxabsdiff_vec(xs1[i], xs2[i])
        if di > dx
            dx = di
            ix = i
        end
    end

    du = 0.0
    iu = 1
    for i in 1:nu
        di = _maxabsdiff_vec(us1[i], us2[i])
        if di > du
            du = di
            iu = i
        end
    end
    return (; dx, ix, du, iu)
end

function _print_trajectory_summary(label, xs, us)
    println("\n[$label] trajectory summary")
    println("  length(xs) = ", length(xs))
    println("  length(us) = ", length(us))
    if !isempty(xs)
        println("  x first = ", xs[1])
        println("  x last  = ", xs[end])
    end
    if !isempty(us)
        println("  u first = ", us[1])
        println("  u last  = ", us[end])
    end
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

function _opts_from_lmi_cfg(cfg::PBRef.LMIConfig; smoke_like::Bool)
    ΔW = smoke_like ? IA.IntervalBox(IA.interval(0.0, 0.0), 1) : cfg.ΔW
    return (
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        λ = cfg.λ,
        rayon_terminal = cfg.rayon_terminal,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = ΔW,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
    )
end

function _make_symbolic_builder(av_mod, params)
    return function (prob, candidate, cfg)
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
end

function _build_new_context(
    problem,
    xs,
    us,
    Ts,
    av_mod,
    params,
    lmi_cfg;
    smoke_like::Bool = false,
)
    nx = length(xs[1])
    nu = length(us[1])
    opts = _opts_from_lmi_cfg(lmi_cfg; smoke_like = smoke_like)
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
        source = :diagnostic,
        metadata = (; nx, nu),
    )
    builder = _make_symbolic_builder(av_mod, params)
    ctx = SC.build_symbolic_context(problem, cand, cert_cfg, builder)
    return (; ctx, cert_cfg, cand)
end

function _linearize(ctx, xk, uk)
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

function _compare_linearization_pair(labelA, ctxA, labelB, ctxB, xs, us, i_state)
    xk = xs[i_state]
    uk = us[i_state]
    affA, LA = _linearize(ctxA, xk, uk)
    affB, LB = _linearize(ctxB, xk, uk)

    dA = _maxabsdiff_mat(affA.A, affB.A)
    dB = _maxabsdiff_mat(affA.B, affB.B)
    dC = _maxabsdiff_vec(affA.c, affB.c)
    dD = _maxabsdiff_mat(affA.D, affB.D)
    dL = _maxabsdiff_vec(LA, LB)

    println("  compare @i=$i_state ($labelA vs $labelB)")
    @printf(
        "    max|ΔA| = %.3e | max|ΔB| = %.3e | max|Δc| = %.3e | max|ΔD| = %.3e | max|ΔL| = %.3e\n",
        dA,
        dB,
        dC,
        dD,
        dL
    )
    return (; dA, dB, dC, dD, dL)
end

function _backend(ctx)
    if hasproperty(ctx, :backend)
        return getproperty(ctx, :backend)
    end
    return getproperty(ctx, :sdp_opt)
end

function _terminal_ellipsoid(xs, radius::Float64)
    nx = length(xs[end])
    P = Matrix{Float64}(I, nx, nx) * (1.0 / radius^2)
    return UT.Ellipsoid(P, _vecf(xs[end]))
end

function _run_backward_prefix(ctx, xs, us; nsteps::Int, radius::Float64)
    n_u = length(us)
    n_use = min(nsteps, n_u)
    E_next = _terminal_ellipsoid(xs, radius)
    out = NamedTuple[]
    for j in 1:n_use
        i_state = n_u - j + 1
        xk = xs[i_state]
        uk = us[i_state]
        affineSys, L = _linearize(ctx, xk, uk)
        E_prev, kappa, cost = SY.transition_backward(
            affineSys,
            E_next,
            _vecf(xk),
            _vecf(uk),
            ctx.Uformat,
            ctx.Wformat,
            ctx.S,
            L,
            _backend(ctx);
            maxδx = ctx.maxδx,
            maxδu = ctx.maxδu,
            λ = ctx.λ,
        )
        status = (E_prev !== nothing && kappa !== nothing) ? :ok : :fail
        push!(out, (; j, i_state, status, cost = cost === nothing ? NaN : Float64(cost)))
        status == :fail && break
        E_next = E_prev
    end
    return out
end

function _print_prefix_result(label, rows)
    println("\n[$label] backward prefix")
    isempty(rows) && (println("  no step"); return)
    for r in rows
        if r.status == :ok
            @printf("  j=%d i=%d status=ok cost=%.6g\n", r.j, r.i_state, r.cost)
        else
            @printf("  j=%d i=%d status=FAIL\n", r.j, r.i_state)
        end
    end
    first_fail = findfirst(x -> x.status == :fail, rows)
    if first_fail === nothing
        println("  prefix result: no failure on tested horizon")
    else
        println(
            "  prefix result: first fail at tested step j=",
            rows[first_fail].j,
            " (state index i=",
            rows[first_fail].i_state,
            ")",
        )
    end
end

function main(args = ARGS)
    parsed = _parse_args(args)
    nsteps = parsed.nsteps
    skip_permis_e = parsed.skip_permis_e

    println("=== LMI mismatch diagnosis ===")
    println("nsteps prefix = ", nsteps)
    println("skip Permis_E nominal = ", skip_permis_e)

    lmi_cfg = PBRef.LMIConfig()
    println("\n[config]")
    println("  λ = ", lmi_cfg.λ)
    println("  maxδx = ", lmi_cfg.maxδx)
    println("  maxδu = ", lmi_cfg.maxδu)
    println("  rayon_terminal = ", lmi_cfg.rayon_terminal)
    println("  ΔW reference (Permis_B) = ", lmi_cfg.ΔW)
    println("  ΔW smoke-like = ", IA.IntervalBox(IA.interval(0.0, 0.0), 1))

    pb_run = PBRef.run_nominal_simulation(;
        use_cache = true,
        force_recompute = false,
        save_plots = false,
        show_animation = false,
    )
    pb_run_lmi = PBRef.prepare_run_result_for_lmi(pb_run; enable_unwrap = true)
    xs_pb = pb_run_lmi.state_list
    us_pb = pb_run_lmi.input_list
    _print_trajectory_summary("Permis_B reference", xs_pb, us_pb)

    pb_problem = _make_problem_from_pb_run(pb_run)
    old_ctx_pb = PBRef.build_symbolic_lmi_context(pb_run_lmi, lmi_cfg)
    new_pb_same = _build_new_context(
        pb_problem,
        xs_pb,
        us_pb,
        pb_run.Δt,
        PBRef.AV,
        pb_run.params,
        lmi_cfg;
        smoke_like = false,
    )
    new_pb_smoke = _build_new_context(
        pb_problem,
        xs_pb,
        us_pb,
        pb_run.Δt,
        PBRef.AV,
        pb_run.params,
        lmi_cfg;
        smoke_like = true,
    )

    println("\n[linearization consistency on SAME trajectory (Permis_B)]")
    isamp = unique(
        clamp.(
            Int.([1, floor(length(us_pb) / 2), length(us_pb)]),
            1,
            max(1, length(us_pb)),
        ),
    )
    for i_state in isamp
        _compare_linearization_pair(
            "old Permis_B",
            old_ctx_pb,
            "new certifier (same opts)",
            new_pb_same.ctx,
            xs_pb,
            us_pb,
            i_state,
        )
        _compare_linearization_pair(
            "old Permis_B",
            old_ctx_pb,
            "new certifier (smoke-like ΔW)",
            new_pb_smoke.ctx,
            xs_pb,
            us_pb,
            i_state,
        )
    end

    old_prefix = _run_backward_prefix(
        old_ctx_pb,
        xs_pb,
        us_pb;
        nsteps = nsteps,
        radius = lmi_cfg.rayon_terminal,
    )
    new_same_prefix = _run_backward_prefix(
        new_pb_same.ctx,
        xs_pb,
        us_pb;
        nsteps = nsteps,
        radius = lmi_cfg.rayon_terminal,
    )
    new_smoke_prefix = _run_backward_prefix(
        new_pb_smoke.ctx,
        xs_pb,
        us_pb;
        nsteps = nsteps,
        radius = lmi_cfg.rayon_terminal,
    )
    _print_prefix_result("old Permis_B", old_prefix)
    _print_prefix_result("new certifier (same trajectory, same opts)", new_same_prefix)
    _print_prefix_result("new certifier (same trajectory, smoke-like ΔW)", new_smoke_prefix)

    if !skip_permis_e
        println("\n[Permis_E nominal trajectory check]")
        pe_nom = PERef.get_nominal_traje()
        pe_nom.candidate === nothing && error("Permis_E nominal candidate is nothing.")
        xs_pe_raw = collect(ST.enum_elems(pe_nom.candidate.x_traj))
        us_pe = collect(ST.enum_elems(pe_nom.candidate.u_traj))
        xs_pe =
            PBRef.unwrap_periodic_state_list(xs_pe_raw, SVector(3, 4), SVector(2pi, 2pi))
        _print_trajectory_summary("Permis_E nominal (unwrapped)", xs_pe, us_pe)

        td = _trajectory_delta(xs_pb, us_pb, xs_pe, us_pe)
        println("\n[trajectory delta: Permis_B vs Permis_E]")
        @printf("  max|Δx| = %.6g at i=%d\n", td.dx, td.ix)
        @printf("  max|Δu| = %.6g at i=%d\n", td.du, td.iu)
        println(
            "  same lengths? x=",
            length(xs_pb) == length(xs_pe),
            " u=",
            length(us_pb) == length(us_pe),
        )
        if 1 <= td.ix <= min(length(xs_pb), length(xs_pe))
            println("  x_pb[i] = ", xs_pb[td.ix])
            println("  x_pe[i] = ", xs_pe[td.ix])
        end
        if 1 <= td.iu <= min(length(us_pb), length(us_pe))
            println("  u_pb[i] = ", us_pb[td.iu])
            println("  u_pe[i] = ", us_pe[td.iu])
        end

        pe_problem = pe_nom.problem
        pe_params = PERef.AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
        pe_params_pb = PBRef.AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
        run_pe = (;
            concrete_system = pe_problem.system,
            params = pe_params_pb,
            Δt = pe_nom.candidate.Ts,
            state_list = xs_pe,
            input_list = us_pe,
        )
        old_ctx_pe = PBRef.build_symbolic_lmi_context(run_pe, lmi_cfg)

        new_pe_same = _build_new_context(
            pe_problem,
            xs_pe,
            us_pe,
            pe_nom.candidate.Ts,
            PERef.AV,
            pe_params,
            lmi_cfg;
            smoke_like = false,
        )
        new_pe_smoke = _build_new_context(
            pe_problem,
            xs_pe,
            us_pe,
            pe_nom.candidate.Ts,
            PERef.AV,
            pe_params,
            lmi_cfg;
            smoke_like = true,
        )

        old_pe_prefix = _run_backward_prefix(
            old_ctx_pe,
            xs_pe,
            us_pe;
            nsteps = nsteps,
            radius = lmi_cfg.rayon_terminal,
        )
        new_pe_same_prefix = _run_backward_prefix(
            new_pe_same.ctx,
            xs_pe,
            us_pe;
            nsteps = nsteps,
            radius = lmi_cfg.rayon_terminal,
        )
        new_pe_smoke_prefix = _run_backward_prefix(
            new_pe_smoke.ctx,
            xs_pe,
            us_pe;
            nsteps = nsteps,
            radius = lmi_cfg.rayon_terminal,
        )
        _print_prefix_result("old helper (Permis_E nominal, same opts)", old_pe_prefix)
        _print_prefix_result(
            "new certifier (Permis_E nominal, same opts)",
            new_pe_same_prefix,
        )
        _print_prefix_result(
            "new certifier (Permis_E nominal, smoke-like opts)",
            new_pe_smoke_prefix,
        )
    end

    return println("\n=== diagnosis end ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
