using LinearAlgebra
using Printf
using JuMP
import IntervalArithmetic as IA
import MathOptInterface as MOI
import StaticArrays: SVector
import MosekTools
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PR = DI.Problem
const ST = DI.System
const OP = DI.Optim
const SC = OP.SymbolicCertifier
const SY = DI.Symbolic

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

Base.@kwdef struct SmokeConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 700
    terminal_radius::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
    )
    ΔU::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-1.0, 1.0), IA.interval(-1.2, 1.2))
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)
    debug_k_start::Union{Nothing, Int} = nothing
    debug_k_stop::Int = 0
end

function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims and periodic_periods must have same length.")
    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        1 <= d <= nx || error("Invalid periodic dimension: $d")
        p > 0 || error("Invalid period <= 0: $p")
        for k in 2:length(xs)
            Δ_raw = xs[k][d] - xs[k - 1][d]
            Δ_min = Δ_raw - round(Δ_raw / p) * p
            xs[k][d] = xs[k - 1][d] + Δ_min
        end
    end
    return [SVector{nx, Float64}(x) for x in xs]
end

function build_concrete_system()
    x_domain = UT.HyperRectangle(SVector(-1.0, -1.0, -pi, -pi), SVector(10.0, 9.0, pi, pi))
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(50.0))
    obstacles_xy = [
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 6.0), SVector(10.0, 9.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)
    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))
    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)
    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_domain()
    inputs_delta = [
        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return DO.CustomList(inputs)
end

function build_control_problem(system)
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    initial_set =
        UT.HyperRectangle(SVector(-1.0, -1.0, -0.4, -0.4), SVector(1.0, 1.0, 0.4, 0.4))
    target_set = UT.HyperRectangle(
        SVector(9.0, 5.0, pi - 5 * (pi / 180), -5 * (pi / 180)),
        SVector(10.0, 6.0, pi + 5 * (pi / 180), 5 * (pi / 180)),
    )
    prob = PR.OptimalControlProblem(
        system,
        initial_set,
        target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
    return (; x0, initial_set, target_set, problem = prob)
end

periodicity_kwargs(cfg::SmokeConfig) = (;
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

function get_nominal_candidate(cfg::SmokeConfig)
    sys = build_concrete_system()
    ctl = build_control_problem(sys.concrete_system)
    Udom = build_input_domain()

    gen_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        Udom,
        AV.jacobian_bound(sys.params),
        periodicity_kwargs(cfg),
        cfg.nstep,
        _ -> ctl.x0,
    )
    gen = OP.CenteredAbstractionGenerator{typeof(ctl.problem), typeof(gen_cfg), Any, Any}(
        ctl.problem,
        gen_cfg,
        nothing,
        nothing,
        false,
        0.0,
    )
    OP.set_problem!(gen, ctl.problem)
    OP.generate!(gen)
    cand = OP.get_trajectory(gen)
    cand === nothing && error("Nominal trajectory generation failed.")
    return (;
        candidate = cand,
        problem = ctl.problem,
        params = sys.params,
        success = OP.get_success(gen),
        solve_time = OP.get_solve_time(gen),
    )
end

function make_backend()
    return optimizer_with_attributes(
        MosekTools.Optimizer,
        MOI.Silent() => true,
        MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => 0,
        MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
    )
end

function make_symbolic_builder(params)
    return function (prob, candidate, cert_cfg)
        o = cert_cfg.options
        return AV.symbolic_system(
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

function ellipsoid_metrics(E)
    P = Matrix{Float64}(E.P)
    evals = eigvals(Symmetric(P))
    λmin = minimum(evals)
    λmax = maximum(evals)
    condP = λmax / λmin
    vol = UT.get_volume(E)
    return (; vol, λmin, λmax, condP)
end

function print_storage_audit(res)
    @assert res !== nothing
    @assert res.lmi_data !== nothing
    ells = res.lmi_data.ellipsoids
    steps = res.steps
    n_ok = count(s -> s.status == :ok, steps)
    n_fail = count(s -> s.status != :ok, steps)

    println("\n=== Storage Audit ===")
    println(
        "n_steps=",
        length(steps),
        " | ok=",
        n_ok,
        " | fail=",
        n_fail,
        " | n_ellipsoids=",
        length(ells),
    )
    println("expect n_ellipsoids = n_ok + 1 -> ", length(ells) == n_ok + 1)
    if n_fail > 0
        first_fail = findfirst(s -> s.status != :ok, steps)
        println("first_failed_step_idx=", first_fail, " | failed_k=", steps[first_fail].k)
    end
    if !isempty(ells)
        mT = ellipsoid_metrics(ells[1])
        @printf(
            "terminal | vol=%.12e | cond(P)=%.3e | λmin=%.3e | λmax=%.3e\n",
            mT.vol,
            mT.condP,
            mT.λmin,
            mT.λmax
        )
    end
    for i in 2:length(ells)
        k = steps[i - 1].k
        m = ellipsoid_metrics(ells[i])
        @printf(
            "k=%d | vol=%.12e | cond(P)=%.3e | λmin=%.3e | λmax=%.3e\n",
            k,
            m.vol,
            m.condP,
            m.λmin,
            m.λmax
        )
    end
end

function run_manual_backward_debug(
    ctx;
    k_start::Union{Nothing, Int} = nothing,
    k_stop::Int = 0,
)
    nx = length(ctx.xs[end])
    PN = Matrix{Float64}(I, nx, nx) * (1.0 / ctx.terminal_radius^2)
    E_next = UT.Ellipsoid(PN, ctx.xs[end])
    K = ctx.K
    ks = k_start === nothing ? (K - 1) : Int(k_start)
    ks = clamp(ks, 0, K - 1)
    ke = clamp(k_stop, 0, ks)

    println("\n=== Manual Backward Debug ===")
    println("K=", K, " | loop k=", ks, ":-1:", ke, " (physical old-index = k+1)")
    mN = ellipsoid_metrics(E_next)
    @printf("start E_next (terminal) | vol=%.12e | cond(P)=%.3e\n", mN.vol, mN.condP)

    for k in ks:-1:ke
        xk = ctx.xs[k + 1]
        uk = ctx.us[k + 1]
        wk = zeros(length(ctx.w))
        Xbar = IA.IntervalBox(xk .+ ctx.ΔX)
        Ubar = IA.IntervalBox(uk .+ ctx.ΔU)
        Wbar = IA.IntervalBox(wk .+ ctx.ΔW)
        affineSys, L = ST.buildAffineApproximation(
            ctx.fsymbolic,
            ctx.x,
            ctx.u,
            ctx.w,
            xk,
            uk,
            wk,
            Xbar,
            Ubar,
            Wbar,
        )
        E_prev, _, cost = SY.transition_backward(
            affineSys,
            E_next,
            xk,
            uk,
            ctx.Uformat,
            ctx.Wformat,
            ctx.S,
            L,
            ctx.backend;
            maxδx = ctx.maxδx,
            maxδu = ctx.maxδu,
            λ = ctx.λ,
        )
        if E_prev === nothing
            println("k=", k, " | status=FAIL")
            break
        end
        m = ellipsoid_metrics(E_prev)
        @printf(
            "k=%d (old=%d) | ū=%s | Lips4=[%.6g, %.6g, %.6g, %.6g] | cost=%.6g | vol=%.12e | cond(P)=%.3e\n",
            k,
            k + 1,
            string(uk),
            L[1],
            L[2],
            L[3],
            L[4],
            Float64(cost),
            m.vol,
            m.condP,
        )
        E_next = E_prev
    end
end

function main(cfg::SmokeConfig = SmokeConfig())
    nominal = get_nominal_candidate(cfg)
    println("generator_success=", nominal.success, " | solve_time=", nominal.solve_time)

    cand_raw = nominal.candidate
    xs_raw = collect(ST.enum_elems(cand_raw.x_traj))
    us = collect(ST.enum_elems(cand_raw.u_traj))
    xs = unwrap_periodic_state_list(xs_raw, cfg.periodic_dims, cfg.periodic_periods)

    cand = OP.CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = cand_raw.Ts,
        source = cand_raw.source,
        metadata = cand_raw.metadata,
    )

    opts = (
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        λ = cfg.λ,
        rayon_terminal = cfg.terminal_radius,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
    )
    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        nominal.problem.system.X,
        nominal.problem.system.U,
        Wdom,
        make_backend(),
        opts,
    )
    builder = make_symbolic_builder(nominal.params)

    cert = SC.EllipsoidalBackwardCertifier{
        typeof(nominal.problem),
        typeof(cand),
        typeof(cert_cfg),
        Any,
        typeof(builder),
    }(
        nothing,
        nothing,
        cert_cfg,
        nothing,
        false,
        0.0,
        builder,
    )

    SC.set_problem!(cert, nominal.problem)
    SC.set_trajectory!(cert, cand)
    SC.certify!(cert)
    res = SC.get_result(cert)
    @assert res !== nothing

    println("certifier_success=", SC.get_success(cert), " | n_steps=", length(res.steps))
    print_storage_audit(res)

    ctx = SC.build_symbolic_context(nominal.problem, cand, cert_cfg, builder)
    return run_manual_backward_debug(
        ctx;
        k_start = cfg.debug_k_start,
        k_stop = cfg.debug_k_stop,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
