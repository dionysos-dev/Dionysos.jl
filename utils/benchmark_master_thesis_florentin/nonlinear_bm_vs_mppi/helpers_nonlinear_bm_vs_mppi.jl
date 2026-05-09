using FillArrays
import GLPK
import HybridSystems
import IntervalArithmetic as IA
using JuMP: optimizer_with_attributes
import LinearAlgebra as LA
import MathOptInterface as MOI
using MathematicalSystems
using Plots
import Polyhedra
import Random
import Statistics
using SemialgebraicSets
import StaticArrays: SVector
using Printf: @sprintf

const HS = HybridSystems
const POLY_LIB = Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer)

Base.@kwdef struct NonlinearBMvsMPPIConfig
    horizons::Tuple{Vararg{Int}} = (6, 8, 10, 12, 14, 16, 18, 20)
    Ts::Float64 = 1.0
    μ::Float64 = 0.00005
    noise::Bool = false

    output_root::String = joinpath(@__DIR__, "outputs")

    eval_step_cost::Float64 = 1.0
    eval_control_weight::Float64 = 0.05
    eval_terminal_weight::Float64 = 25.0
    eval_terminal_miss_penalty::Float64 = 1.0e4
    eval_bad_cost::Float64 = 1.0e8

    bm_target_inner_scale::Float64 = 0.95
    bm_state_target_weight::Float64 = 0.05
    bm_indicator::Bool = false

    mppi_nsamples::Int = 1000
    mppi_niter::Int = 12
    mppi_λ::Float64 = 5.0
    mppi_noise_u::Float64 = 3.0
    mppi_greedy_gain::Float64 = 0.85
    mppi_warm_start_mode::Symbol = :greedy_target

    λ::Float64 = 0.05
    terminal_radius::Float64 = 1.0
    maxδx::Float64 = 100.0
    maxδu::Float64 = 60.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-2.0, 2.0),
        IA.interval(-2.0, 2.0),
    )
    ΔW::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(0.0, 0.0),
        IA.interval(0.0, 0.0),
    )

    rng_seed::Int = 1234
    verbose::Bool = false
end

Base.@kwdef struct NeutralSeedConfig
    nstep::Int
    warm_start_mode::Symbol = :greedy_target
    greedy_gain::Float64 = 0.85
end

mutable struct NeutralSeedGenerator{P, C} <: OP.AbstractHeuristicGenerator
    problem::Union{Nothing, P}
    config::C
    candidate::Union{Nothing, OP.CandidateTrajectory}
    success::Bool
    solve_time_sec::Float64
end

function NeutralSeedGenerator(problem, config::NeutralSeedConfig)
    return NeutralSeedGenerator{Any, typeof(config)}(problem, config, nothing, false, 0.0)
end

function OP.set_problem!(gen::NeutralSeedGenerator, problem)
    gen.problem = problem
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    return gen
end

OP.get_trajectory(gen::NeutralSeedGenerator) = gen.candidate
OP.get_success(gen::NeutralSeedGenerator) = gen.success
OP.get_solve_time(gen::NeutralSeedGenerator) = gen.solve_time_sec

function OP.generate!(gen::NeutralSeedGenerator)
    gen.problem === nothing && error("Call set_problem!(seed_generator, problem) first.")
    t0 = time()
    controls = build_seed_controls(gen.problem, gen.config)
    gen.candidate = candidate_from_controls(
        gen.problem,
        controls;
        source = :neutral_seed,
        metadata = (; warm_start_mode = gen.config.warm_start_mode),
    )
    gen.success = final_state(gen.candidate) ∈ gen.problem.target_set
    gen.solve_time_sec = time() - t0
    return gen
end

function benchmark_output_paths(cfg::NonlinearBMvsMPPIConfig)
    root = cfg.output_root
    data_dir = joinpath(root, "data")
    plots_dir = joinpath(root, "plots")
    summary_dir = joinpath(root, "summary")
    mkpath(data_dir)
    mkpath(plots_dir)
    mkpath(summary_dir)
    return (; root, data_dir, plots_dir, summary_dir)
end

function horizon_paths(cfg::NonlinearBMvsMPPIConfig, N::Int)
    base = joinpath(cfg.output_root, @sprintf("horizon_%02d", N))
    data_dir = joinpath(base, "data")
    plots_dir = joinpath(base, "plots")
    mkpath(data_dir)
    mkpath(plots_dir)
    return (; root = base, data_dir, plots_dir)
end

function build_nonlinear_finite_horizon_problem(
    cfg::NonlinearBMvsMPPIConfig,
    N::Int;
    nonlinear_module = NL,
)
    base_problem = nonlinear_module.problem(; Ts = cfg.Ts, μ = cfg.μ, noise = cfg.noise)
    return PR.OptimalControlProblem(
        base_problem.system,
        base_problem.initial_set,
        base_problem.target_set,
        base_problem.state_cost,
        base_problem.transition_cost,
        N,
    )
end

function domain_rectangle(problem)
    X = problem.system.X
    return UT.HyperRectangle(
        SVector(IA.inf(X[1]), IA.inf(X[2])),
        SVector(IA.sup(X[1]), IA.sup(X[2])),
    )
end

function project_input_to_domain(u, Uset)
    uvec = SVector{length(u), Float64}(u)
    if Uset isa UT.HyperRectangle
        return SVector{length(uvec), Float64}(clamp.(uvec, Uset.lb, Uset.ub))
    end
    return uvec
end

zero_noise(problem) = zeros(Float64, UT.get_dims(problem.system.W))
initial_state(problem) = SVector{2, Float64}(UT.get_center(problem.initial_set))
target_center(problem) = SVector{2, Float64}(UT.get_center(problem.target_set))

function discrete_rollout_step(problem, x, u)
    noise = zero_noise(problem)
    return SVector{2, Float64}(problem.system.f_eval(collect(x), collect(u), noise))
end

final_state(cand::OP.CandidateTrajectory) = collect(ST.enum_elems(cand.x_traj))[end]

function target_distance(target_set::UT.Ellipsoid, x)
    dx = SVector{length(x), Float64}(x) - target_set.c
    value = sqrt(max(0.0, LA.dot(dx, target_set.P * dx)))
    return max(0.0, value - 1.0)
end

function obstacle_hit(problem, x)
    if !hasproperty(problem.system, :obstacles)
        return false
    end
    return any(obs -> (x ∈ obs), problem.system.obstacles)
end

function common_rollout_cost(problem, xs, us, cfg::NonlinearBMvsMPPIConfig)
    isempty(xs) && return cfg.eval_bad_cost

    total_cost = 0.0

    for k in eachindex(us)
        xk = xs[k]
        uk = us[k]
        if !(xk ∈ problem.system.X) || obstacle_hit(problem, xk)
            return total_cost + cfg.eval_bad_cost
        end
        total_cost += cfg.eval_step_cost
        total_cost += cfg.eval_control_weight * sum(abs2, uk)
    end

    xT = xs[end]
    if !(xT ∈ problem.system.X) || obstacle_hit(problem, xT)
        return total_cost + cfg.eval_bad_cost
    end

    distT = target_distance(problem.target_set, xT)
    total_cost += cfg.eval_terminal_weight * distT^2
    if !(xT ∈ problem.target_set)
        total_cost += cfg.eval_terminal_miss_penalty
    end

    return total_cost
end

function rollout_controls(problem, controls; truncate_on_target::Bool)
    x = initial_state(problem)
    xs = SVector{2, Float64}[x]
    us = SVector{2, Float64}[]

    for u in controls
        uproj = project_input_to_domain(u, problem.system.U)
        push!(us, uproj)
        x = discrete_rollout_step(problem, x, uproj)
        push!(xs, x)
        if truncate_on_target && (x ∈ problem.target_set)
            break
        end
    end

    return xs, us
end

function candidate_from_controls(problem, controls; source::Symbol, metadata = (;))
    xs, us = rollout_controls(problem, controls; truncate_on_target = true)
    return OP.CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = problem.system.Ts,
        source = source,
        metadata = metadata,
    )
end

function evaluate_true_rollout(problem, candidate, cfg::NonlinearBMvsMPPIConfig)
    if candidate === nothing
        return (
            success = false,
            first_hit_step = nothing,
            realized_horizon = 0,
            true_rollout_cost = NaN,
            domain_violation = true,
            obstacle_violation = false,
            terminal_distance = NaN,
        )
    end

    controls = collect(ST.enum_elems(candidate.u_traj))
    xs, us = rollout_controls(problem, controls; truncate_on_target = false)
    hit_idx = findfirst(x -> (x ∈ problem.target_set), xs)
    xT = xs[end]

    return (
        success = hit_idx !== nothing,
        first_hit_step = isnothing(hit_idx) ? nothing : hit_idx - 1,
        realized_horizon = length(us),
        true_rollout_cost = common_rollout_cost(problem, xs, us, cfg),
        domain_violation = any(x -> !(x ∈ problem.system.X), xs),
        obstacle_violation = any(x -> obstacle_hit(problem, x), xs),
        terminal_distance = target_distance(problem.target_set, xT),
    )
end

function greedy_target_control(problem, x, gain::Float64)
    xvec = SVector{2, Float64}(x)
    drift = discrete_rollout_step(problem, xvec, SVector(0.0, 0.0))
    u = gain .* ((target_center(problem) .- drift) ./ problem.system.Ts)
    return project_input_to_domain(u, problem.system.U)
end

function build_seed_controls(problem, seed_cfg::NeutralSeedConfig)
    controls = SVector{2, Float64}[]
    x = initial_state(problem)

    for _ in 1:seed_cfg.nstep
        u = if seed_cfg.warm_start_mode == :zero
            SVector(0.0, 0.0)
        elseif seed_cfg.warm_start_mode == :greedy_target
            greedy_target_control(problem, x, seed_cfg.greedy_gain)
        else
            error("Unsupported warm_start_mode=$(seed_cfg.warm_start_mode)")
        end
        push!(controls, u)
        x = discrete_rollout_step(problem, x, u)
        x ∈ problem.target_set && break
    end

    if isempty(controls)
        push!(controls, SVector(0.0, 0.0))
    end
    return controls
end

function build_nonlinear_mppi_generator(problem, cfg::NonlinearBMvsMPPIConfig, N::Int)
    seed_cfg = NeutralSeedConfig(
        nstep = N,
        warm_start_mode = cfg.mppi_warm_start_mode,
        greedy_gain = cfg.mppi_greedy_gain,
    )
    seed_gen = NeutralSeedGenerator(problem, seed_cfg)

    discrete_dynamics = (prob, x, u, _k, _Δt) -> discrete_rollout_step(prob, x, u)
    noise_sampler = (rng, _u, _k) -> SVector(
        cfg.mppi_noise_u * Random.randn(rng),
        cfg.mppi_noise_u * Random.randn(rng),
    )
    project_input = u -> project_input_to_domain(u, problem.system.U)

    trajectory_cost =
        (prob, cand) -> common_rollout_cost(
            prob,
            collect(ST.enum_elems(cand.x_traj)),
            collect(ST.enum_elems(cand.u_traj)),
            cfg,
        )
    success_fun = (prob, cand) -> final_state(cand) ∈ prob.target_set

    mppi_cfg = OP.MPPIConfig(
        problem.system.Ts,
        N,
        cfg.mppi_nsamples,
        cfg.mppi_niter,
        cfg.mppi_λ,
        initial_state,
        discrete_dynamics,
        noise_sampler,
        project_input,
        trajectory_cost,
        success_fun,
        (prob, x) -> SVector{2, Float64}(x),
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function solve_mppi_candidate(problem, cfg::NonlinearBMvsMPPIConfig, N::Int)
    gen = build_nonlinear_mppi_generator(problem, cfg, N)
    OP.set_problem!(gen, problem)
    Random.seed!(cfg.rng_seed + 1000 * N)
    OP.generate!(gen)

    return (
        method = :mppi,
        candidate = OP.get_trajectory(gen),
        seed_candidate = OP.get_seed(gen),
        success = OP.get_success(gen),
        solve_time_sec = OP.get_solve_time(gen),
        diagnostics = gen.diagnostics,
        solver_status = OP.get_success(gen) ? "generated" : "candidate_generated",
        surrogate_used = false,
        surrogate_note = "MPPI uses the true nonlinear discrete rollout, with a neutral deterministic warm start and no abstraction reference tracking.",
    )
end

function inscribed_target_box(problem, cfg::NonlinearBMvsMPPIConfig)
    radii = UT.get_length_semiaxis(problem.target_set)
    r = cfg.bm_target_inner_scale * minimum(radii) / sqrt(2.0)
    c = target_center(problem)
    return UT.HyperRectangle(c .- r, c .+ r)
end

function obstacle_outer_box(problem)
    obs = first(problem.system.obstacles)
    box = UT.get_min_bounding_box(obs)
    return UT.HyperRectangle(
        SVector(IA.inf(box[1]), IA.inf(box[2])),
        SVector(IA.sup(box[1]), IA.sup(box[2])),
    )
end

function has_positive_area_intersection(a::UT.HyperRectangle, b::UT.HyperRectangle; tol = 1.0e-9)
    return all(max.(a.lb, b.lb) .< (min.(a.ub, b.ub) .- tol))
end

function coarse_surrogate_cells(problem, obs_box)
    dom = domain_rectangle(problem)
    xL, xU = dom.lb[1], dom.ub[1]
    yL, yU = dom.lb[2], dom.ub[2]
    xoL, xoU = obs_box.lb[1], obs_box.ub[1]
    yoL, yoU = obs_box.lb[2], obs_box.ub[2]

    cells = [
        UT.HyperRectangle(SVector(xL, yL), SVector(xoL, yoL)),
        UT.HyperRectangle(SVector(xL, yoL), SVector(xoL, yU)),
        UT.HyperRectangle(SVector(xoL, yL), SVector(xU, yoL)),
        UT.HyperRectangle(SVector(xoL, yoU), SVector(xoU, yU)),
        UT.HyperRectangle(SVector(xoU, yoL), SVector(xU, yoU)),
        UT.HyperRectangle(SVector(xoU, yoU), SVector(xU, yU)),
    ]

    return [rect for rect in cells if !isempty(rect) && !has_positive_area_intersection(rect, obs_box)]
end

function polyhedron_from_rectangle(rect::UT.HyperRectangle)
    rep = (
        Polyhedra.HalfSpace(SVector(1.0, 0.0), rect.ub[1]) ∩
        Polyhedra.HalfSpace(SVector(-1.0, 0.0), -rect.lb[1]) ∩
        Polyhedra.HalfSpace(SVector(0.0, 1.0), rect.ub[2]) ∩
        Polyhedra.HalfSpace(SVector(0.0, -1.0), -rect.lb[2])
    )
    return Polyhedra.polyhedron(rep, POLY_LIB)
end

function shared_linear_map(cfg::NonlinearBMvsMPPIConfig)
    A = [
        1.1 -0.2
        0.2 1.1
    ]
    B = cfg.Ts .* Matrix{Float64}(LA.I, 2, 2)
    return A, B
end

function affine_image_box(
    rect::UT.HyperRectangle,
    Uset::UT.HyperRectangle,
    A::AbstractMatrix,
    B::AbstractMatrix,
)
    points = SVector{2, Float64}[]
    for xv in UT.collect_vertices(rect), uv in UT.collect_vertices(Uset)
        push!(points, SVector{2, Float64}(A * xv + B * uv))
    end
    lb = SVector(minimum(p[1] for p in points), minimum(p[2] for p in points))
    ub = SVector(maximum(p[1] for p in points), maximum(p[2] for p in points))
    return UT.HyperRectangle(lb, ub)
end

function build_bm_surrogate_problem(problem, cfg::NonlinearBMvsMPPIConfig, N::Int)
    x0 = initial_state(problem)
    Uset = problem.system.U
    target_box = inscribed_target_box(problem, cfg)
    obs_box = obstacle_outer_box(problem)
    free_cells = coarse_surrogate_cells(problem, obs_box)

    q0 = findfirst(rect -> (x0 ∈ rect), free_cells)
    q0 === nothing && error("Initial point $(x0) is not covered by the BM surrogate partition.")

    target_mode = length(free_cells) + 1
    mode_rects = vcat(free_cells, [target_box])
    mode_polys = [polyhedron_from_rectangle(rect) for rect in mode_rects]

    automaton = HS.GraphAutomaton(length(mode_rects))
    transition_pairs = Tuple{Int, Int}[]

    input_poly = polyhedron_from_rectangle(Uset)
    A_shared, B_shared = shared_linear_map(cfg)
    for i in eachindex(free_cells)
        reach_box = affine_image_box(free_cells[i], Uset, A_shared, B_shared)

        for j in eachindex(mode_rects)
            UT.is_intersection(reach_box, mode_rects[j]) || continue
            HS.add_transition!(automaton, i, j, 1)
            push!(transition_pairs, (i, j))
        end
    end

    isempty(transition_pairs) && error("The BM surrogate has no transitions.")

    resetmaps = [
        MathematicalSystems.ConstrainedLinearControlMap(
            A_shared,
            B_shared,
            FullSpace(),
            input_poly,
        ),
    ]
    switchings = [HS.ControlledSwitching()]

    hybrid_system = HS.HybridSystem(
        automaton,
        [MathematicalSystems.ConstrainedContinuousIdentitySystem(2, p) for p in mode_polys],
        resetmaps,
        switchings,
    )

    mode_costs = [
        UT.ConstantFunction(
            cfg.eval_step_cost +
            cfg.bm_state_target_weight * sum(abs2, UT.get_center(rect) - target_center(problem)),
        ) for rect in free_cells
    ]
    push!(mode_costs, UT.ConstantFunction(0.0))
    trans_costs = [UT.ZeroFunction()]
    modes_schedule = [collect(1:length(free_cells)) for _ in 1:max(N - 1, 0)]
    push!(modes_schedule, [target_mode])

    bm_problem = PR.OptimalControlProblem(
        hybrid_system,
        (q0, collect(x0)),
        target_mode,
        Fill(mode_costs, N),
        Fill(trans_costs, N),
        N,
    )

    metadata = (
        target_box = target_box,
        obstacle_outer_box = obs_box,
        free_cells = free_cells,
        q0 = q0,
        target_mode = target_mode,
        transition_pairs = transition_pairs,
        modes_schedule = modes_schedule,
    )
    return bm_problem, metadata
end

function candidate_from_bm_trajectory(problem, traj; metadata = (;))
    controls = [SVector{2, Float64}(u) for u in traj.u]
    return candidate_from_controls(problem, controls; source = :bm, metadata = metadata)
end

function solve_bm_candidate(problem, cfg::NonlinearBMvsMPPIConfig, N::Int)
    bm_problem, surrogate = build_bm_surrogate_problem(problem, cfg, N)

    milp_solver = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true)
    optimizer_factory = optimizer_with_attributes(
        OP.BemporadMorari.Optimizer{Float64},
        "continuous_solver" => milp_solver,
        "mixed_integer_solver" => milp_solver,
        "indicator" => cfg.bm_indicator,
        "log_level" => 0,
        "problem" => bm_problem,
        "modes" => surrogate.modes_schedule,
    )

    bm_optimizer = MOI.instantiate(optimizer_factory)
    MOI.optimize!(bm_optimizer)

    term = MOI.get(bm_optimizer, MOI.TerminationStatus())
    primal = MOI.get(bm_optimizer, MOI.PrimalStatus())
    raw_status = try
        MOI.get(bm_optimizer, MOI.RawStatusString())
    catch
        ""
    end

    candidate = nothing
    if term == MOI.OPTIMAL || primal == MOI.FEASIBLE_POINT
        traj = MOI.get(bm_optimizer, ST.ContinuousTrajectoryAttribute())
        candidate = candidate_from_bm_trajectory(
            problem,
            traj;
            metadata = (
                ;
                surrogate_transition_count = length(surrogate.transition_pairs),
                surrogate_mode_count = length(surrogate.free_cells) + 1,
            ),
        )
    end

    return (
        method = :bm,
        candidate = candidate,
        seed_candidate = nothing,
        success = candidate !== nothing,
        solve_time_sec = MOI.get(bm_optimizer, MOI.SolveTimeSec()),
        diagnostics = (
            termination_status = string(term),
            primal_status = string(primal),
            raw_status = raw_status,
        ),
        solver_status = string(term),
        surrogate_used = true,
        surrogate_note = "BM is solved as a MILP on a benchmark-local PWA/hybrid surrogate with local linear maps, a conservative obstacle box, an inner target box, and mode-dependent constant stage costs; the resulting control sequence is always re-simulated on the true nonlinear dynamics before evaluation and certification.",
        surrogate = surrogate,
    )
end

function build_nonlinear_symbolic_builder()
    return function (prob, _candidate, certifier_cfg)
        o = certifier_cfg.options
        sys = prob.system
        return ST.SymbolicSystem(
            sys.fsymbolicT,
            sys.fsymbolic,
            sys.Ts,
            sys.nx,
            sys.nu,
            sys.nw,
            sys.x,
            sys.u,
            sys.w,
            o.ΔX,
            o.ΔU,
            o.ΔW,
            sys.X,
            sys.U,
            sys.W,
            sys.obstacles,
            sys.f_eval,
            sys.f_backward_eval,
            sys.Uformat,
            sys.Wformat,
        )
    end
end

function build_nonlinear_certifier(problem, cfg::NonlinearBMvsMPPIConfig)
    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
    )
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        problem.system.W,
        build_backend(; verbose = cfg.verbose),
        opts,
    )
    return SC.EllipsoidalBackwardCertifier(cert_cfg, build_nonlinear_symbolic_builder())
end

function compute_ellipsoid_area_metrics(ellipsoids)
    areas = [UT.get_volume(E) for E in ellipsoids]
    if isempty(areas)
        return (
            count = 0,
            mean = NaN,
            median = NaN,
            min = NaN,
            max = NaN,
            sum = NaN,
            areas = Float64[],
        )
    end
    return (
        count = length(areas),
        mean = Statistics.mean(areas),
        median = Statistics.median(areas),
        min = minimum(areas),
        max = maximum(areas),
        sum = sum(areas),
        areas = areas,
    )
end

function certify_candidate(problem, candidate, cfg::NonlinearBMvsMPPIConfig)
    if candidate === nothing
        return (
            attempted = false,
            success = false,
            solve_time_sec = NaN,
            certified_steps = 0,
            full_certified = false,
            failed_k = nothing,
            ellipsoids = UT.Ellipsoid[],
            area_metrics = compute_ellipsoid_area_metrics(UT.Ellipsoid[]),
            result = nothing,
            error_message = "",
        )
    end

    cert = build_nonlinear_certifier(problem, cfg)
    try
        SC.set_problem!(cert, problem)
        SC.set_trajectory!(cert, candidate)
        SC.certify!(cert)
        result = SC.get_result(cert)
        ellipsoids = extract_ellipsoids(result; max_keep = 10_000)
        certified_steps =
            result === nothing ? 0 : count(step -> step.status == :ok, result.steps)
        return (
            attempted = true,
            success = SC.get_success(cert),
            solve_time_sec = SC.get_solve_time(cert),
            certified_steps = certified_steps,
            full_certified = SC.get_success(cert),
            failed_k = result === nothing ? nothing : result.failed_k,
            ellipsoids = ellipsoids,
            area_metrics = compute_ellipsoid_area_metrics(ellipsoids),
            result = result,
            error_message = "",
        )
    catch err
        return (
            attempted = true,
            success = false,
            solve_time_sec = NaN,
            certified_steps = 0,
            full_certified = false,
            failed_k = nothing,
            ellipsoids = UT.Ellipsoid[],
            area_metrics = compute_ellipsoid_area_metrics(UT.Ellipsoid[]),
            result = nothing,
            error_message = sprint(showerror, err),
        )
    end
end

function compute_common_prefix_metrics(area_a::AbstractVector, area_b::AbstractVector)
    common_count = min(length(area_a), length(area_b))
    if common_count == 0
        return (common_count = 0, common_steps = 0, mean = NaN, sum = NaN)
    end
    return (
        common_count = common_count,
        common_steps = max(common_count - 1, 0),
        mean = Statistics.mean(area_a[1:common_count]),
        sum = sum(area_a[1:common_count]),
    )
end

function method_color(method::Symbol)
    return method == :bm ? :red : :blue
end

method_label(method::Symbol) = method == :bm ? "BM" : "MPPI"

function plot_problem_geometry!(fig, problem)
    dom = domain_rectangle(problem)
    plot!(fig, dom; color = :grey80, alpha = 0.1, label = "")
    plot!(fig, problem.initial_set; color = :green, alpha = 0.25, label = "initial")
    plot!(fig, problem.target_set; color = :red, alpha = 0.25, label = "target")
    for (idx, obs) in enumerate(problem.system.obstacles)
        plot!(fig, obs; color = :black, alpha = 0.25, label = idx == 1 ? "obstacle" : "")
    end
    xlims!(fig, dom.lb[1], dom.ub[1])
    ylims!(fig, dom.lb[2], dom.ub[2])
    return fig
end

function trajectory_xy(candidate)
    xs = collect(ST.enum_elems(candidate.x_traj))
    return [x[1] for x in xs], [x[2] for x in xs]
end

function control_series(candidate, i::Int)
    us = collect(ST.enum_elems(candidate.u_traj))
    return [u[i] for u in us]
end

function save_empty_plot(path::AbstractString, title::AbstractString, message::AbstractString)
    fig = plot(; legend = false, title = title)
    annotate!(fig, 0.5, 0.5, text(message, 12, :black, :center))
    savefig(fig, path)
    return path
end

function save_phase_comparison_plot!(path, problem, bm_result, mppi_result)
    fig = plot(; aspect_ratio = :equal, legend = :bottomright, title = "Nonlinear benchmark")
    plot_problem_geometry!(fig, problem)

    if mppi_result.seed_candidate !== nothing
        x, y = trajectory_xy(mppi_result.seed_candidate)
        plot!(fig, x, y; lw = 1.5, ls = :dash, color = :steelblue, label = "MPPI seed")
    end

    for result in (bm_result, mppi_result)
        result.candidate === nothing && continue
        x, y = trajectory_xy(result.candidate)
        plot!(
            fig,
            x,
            y;
            lw = 2.5,
            marker = :circle,
            ms = 3,
            color = method_color(result.method),
            label = method_label(result.method),
        )
    end

    savefig(fig, path)
    return path
end

function save_state_control_plot!(path, problem, result)
    result.candidate === nothing &&
        return save_empty_plot(path, "$(method_label(result.method)) states / controls", "No candidate trajectory")

    xs = collect(ST.enum_elems(result.candidate.x_traj))
    us = collect(ST.enum_elems(result.candidate.u_traj))
    tx = (0:(length(xs) - 1)) .* result.candidate.Ts
    tu = (0:(length(us) - 1)) .* result.candidate.Ts
    col = method_color(result.method)

    p1 = plot(tx, [x[1] for x in xs]; color = col, lw = 2, label = "x1", title = "state x1")
    p2 = plot(tx, [x[2] for x in xs]; color = col, lw = 2, label = "x2", title = "state x2")
    u_bounds = problem.system.U
    p3 = plot(
        tu,
        [u[1] for u in us];
        color = col,
        lw = 2,
        label = "u1",
        title = "control u1",
        seriestype = :steppost,
        ylim = (u_bounds.lb[1], u_bounds.ub[1]),
    )
    p4 = plot(
        tu,
        [u[2] for u in us];
        color = col,
        lw = 2,
        label = "u2",
        title = "control u2",
        seriestype = :steppost,
        ylim = (u_bounds.lb[2], u_bounds.ub[2]),
    )

    fig = plot(p1, p2, p3, p4; layout = (2, 2), size = (1000, 700))
    savefig(fig, path)
    return path
end

function save_certification_plot!(path, problem, result, cert_data)
    result.candidate === nothing &&
        return save_empty_plot(path, "$(method_label(result.method)) certification", "No candidate trajectory")

    fig = plot(; aspect_ratio = :equal, legend = :bottomright, title = "$(method_label(result.method)) certification")
    plot_problem_geometry!(fig, problem)
    x, y = trajectory_xy(result.candidate)
    plot!(
        fig,
        x,
        y;
        lw = 2.5,
        marker = :circle,
        ms = 3,
        color = method_color(result.method),
        label = method_label(result.method),
    )

    for (idx, E) in enumerate(cert_data.ellipsoids)
        plot!(
            fig,
            E;
            dims = [1, 2],
            color = :orange,
            alpha = 0.18,
            lw = 1.0,
            label = idx == 1 ? "certified ellipsoids" : "",
        )
    end

    savefig(fig, path)
    return path
end

function row_value(value)
    if value isa AbstractFloat
        return isfinite(value) ? @sprintf("%.6f", value) : ""
    elseif value === nothing
        return ""
    else
        return string(value)
    end
end

function csv_escape(value)
    s = row_value(value)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_rows_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    isempty(rows) && return open(io -> nothing, path, "w")
    cols = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for row in rows
            println(io, join((csv_escape(getproperty(row, c)) for c in cols), ","))
        end
    end
    return path
end

function write_rows_markdown(path::AbstractString, rows::Vector{<:NamedTuple})
    isempty(rows) && return open(io -> nothing, path, "w")
    cols = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, "| ", join(string.(cols), " | "), " |")
        println(io, "| ", join(fill("---", length(cols)), " | "), " |")
        for row in rows
            println(io, "| ", join((replace(row_value(getproperty(row, c)), "|" => "/") for c in cols), " | "), " |")
        end
    end
    return path
end

function method_summary_row(problem, cfg, N, planner_result, eval_data, cert_data, common_metrics)
    diagnostics_string = join(
        ["$(k)=$(v)" for (k, v) in pairs(planner_result.diagnostics)],
        "; ",
    )
    return (
        method = String(method_label(planner_result.method)),
        horizon = N,
        internal_model = planner_result.surrogate_used ? "PWA surrogate for BM" : "true nonlinear discrete rollout",
        solve_time_sec = planner_result.solve_time_sec,
        solver_status = planner_result.solver_status,
        planner_success = planner_result.success,
        true_rollout_success = eval_data.success,
        first_hit_step = eval_data.first_hit_step,
        planned_horizon = N,
        realized_horizon = eval_data.realized_horizon,
        true_rollout_cost = eval_data.true_rollout_cost,
        domain_violation = eval_data.domain_violation,
        obstacle_violation = eval_data.obstacle_violation,
        terminal_distance = eval_data.terminal_distance,
        certification_attempted = cert_data.attempted,
        certification_success = cert_data.success,
        certification_solve_time_sec = cert_data.solve_time_sec,
        certified_steps = cert_data.certified_steps,
        certified_ellipsoid_count = cert_data.area_metrics.count,
        ellipsoid_mean_area = cert_data.area_metrics.mean,
        ellipsoid_median_area = cert_data.area_metrics.median,
        ellipsoid_min_area = cert_data.area_metrics.min,
        ellipsoid_max_area = cert_data.area_metrics.max,
        ellipsoid_sum_area = cert_data.area_metrics.sum,
        common_prefix_steps = common_metrics.common_steps,
        common_prefix_ellipsoid_count = common_metrics.common_count,
        common_prefix_mean_area = common_metrics.mean,
        common_prefix_sum_area = common_metrics.sum,
        cert_failed_k = cert_data.failed_k,
        cert_error = cert_data.error_message,
        note = planner_result.surrogate_note,
        diagnostics = diagnostics_string,
        eval_step_cost = cfg.eval_step_cost,
        eval_control_weight = cfg.eval_control_weight,
        eval_terminal_weight = cfg.eval_terminal_weight,
    )
end

function save_summary_plots!(paths, rows::Vector{<:NamedTuple})
    isempty(rows) && return nothing

    methods = unique(row.method for row in rows)
    horizons = sort(unique(row.horizon for row in rows))

    function metric_series(metric::Symbol, method::String)
        ordered = [row for row in rows if row.method == method]
        sort!(ordered; by = row -> row.horizon)
        xs = [row.horizon for row in ordered]
        ys = [Float64(getproperty(row, metric)) for row in ordered]
        return xs, ys
    end

    fig1 = plot(layout = (2, 2), size = (1100, 800))
    for method in methods
        xs, ys = metric_series(:solve_time_sec, method)
        plot!(fig1[1], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:true_rollout_cost, method)
        plot!(fig1[2], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:true_rollout_success, method)
        plot!(fig1[3], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:certified_steps, method)
        plot!(fig1[4], xs, ys; marker = :circle, lw = 2, label = method)
    end
    plot!(fig1[1]; title = "Solve time", xlabel = "horizon", ylabel = "s")
    plot!(fig1[2]; title = "True rollout cost", xlabel = "horizon")
    plot!(fig1[3]; title = "True rollout success", xlabel = "horizon", ylabel = "0/1", yticks = [0, 1])
    plot!(fig1[4]; title = "Certified steps", xlabel = "horizon")
    savefig(fig1, joinpath(paths.summary_dir, "comparison_core_metrics.pdf"))

    fig2 = plot(layout = (2, 2), size = (1100, 800))
    for method in methods
        xs, ys = metric_series(:ellipsoid_mean_area, method)
        plot!(fig2[1], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:common_prefix_mean_area, method)
        plot!(fig2[2], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:ellipsoid_sum_area, method)
        plot!(fig2[3], xs, ys; marker = :circle, lw = 2, label = method)
        xs, ys = metric_series(:common_prefix_sum_area, method)
        plot!(fig2[4], xs, ys; marker = :circle, lw = 2, label = method)
    end
    plot!(fig2[1]; title = "Mean ellipsoid area", xlabel = "horizon")
    plot!(fig2[2]; title = "Common-prefix mean area", xlabel = "horizon")
    plot!(fig2[3]; title = "Sum ellipsoid area", xlabel = "horizon")
    plot!(fig2[4]; title = "Common-prefix sum area", xlabel = "horizon")
    savefig(fig2, joinpath(paths.summary_dir, "comparison_certification_metrics.pdf"))

    return nothing
end

function save_notes!(paths, cfg::NonlinearBMvsMPPIConfig)
    note_path = joinpath(paths.summary_dir, "benchmark_notes.md")
    lines = [
        "# Nonlinear BM vs MPPI benchmark",
        "",
        "- Dynamics, domain, obstacle, initial ellipsoid, target ellipsoid and input bounds are reused from `problems/non_linear.jl`.",
        "- The finite horizon is overridden locally inside the benchmark code; the original problem file is left untouched.",
        "- Both planners start from the center of the common initial ellipsoid for nominal trajectory generation.",
        "- MPPI uses the true nonlinear discrete rollout and a neutral deterministic warm start, without abstraction or reference tracking.",
        "- BM cannot directly optimize the repository nonlinear symbolic system, so it is run on a benchmark-local PWA/hybrid surrogate only.",
        "- The BM surrogate uses the same input bounds and horizon, local linear maps on a coarse corridor decomposition, a conservative obstacle outer box, an inner target box contained in the true target ellipsoid, and mode-dependent constant stage costs so it remains solvable with the MILP solvers available in this environment.",
        "- Every BM control sequence is re-simulated on the true nonlinear dynamics before common evaluation and certification.",
        "- Both methods are evaluated with the same true-rollout cost and the same ellipsoidal backward certifier configuration.",
        "- Ellipsoid areas are computed with `UT.get_volume`, which matches the repository ellipsoid convention `(x-c)'P(x-c) <= 1`.",
        "- Because the certifier works backward from the terminal state, common-prefix ellipsoid metrics compare the common terminal-nearest certified ellipsoids available for both methods.",
        "",
        "## Default reproducibility settings",
        "",
        "- horizons = $(cfg.horizons)",
        "- rng_seed = $(cfg.rng_seed)",
        "- mppi_nsamples = $(cfg.mppi_nsamples)",
        "- mppi_niter = $(cfg.mppi_niter)",
        "- mppi_lambda = $(cfg.mppi_λ)",
    ]
    open(note_path, "w") do io
        println(io, join(lines, "\n"))
    end
    return note_path
end

function run_single_horizon(cfg::NonlinearBMvsMPPIConfig, N::Int)
    problem = build_nonlinear_finite_horizon_problem(cfg, N)
    paths = horizon_paths(cfg, N)

    bm_result = solve_bm_candidate(problem, cfg, N)
    mppi_result = solve_mppi_candidate(problem, cfg, N)

    bm_eval = evaluate_true_rollout(problem, bm_result.candidate, cfg)
    mppi_eval = evaluate_true_rollout(problem, mppi_result.candidate, cfg)

    bm_cert = certify_candidate(problem, bm_result.candidate, cfg)
    mppi_cert = certify_candidate(problem, mppi_result.candidate, cfg)

    bm_common = compute_common_prefix_metrics(
        bm_cert.area_metrics.areas,
        mppi_cert.area_metrics.areas,
    )
    mppi_common = compute_common_prefix_metrics(
        mppi_cert.area_metrics.areas,
        bm_cert.area_metrics.areas,
    )

    rows = [
        method_summary_row(problem, cfg, N, bm_result, bm_eval, bm_cert, bm_common),
        method_summary_row(problem, cfg, N, mppi_result, mppi_eval, mppi_cert, mppi_common),
    ]

    save_phase_comparison_plot!(
        joinpath(paths.plots_dir, "comparison_phase.pdf"),
        problem,
        bm_result,
        mppi_result,
    )
    save_state_control_plot!(joinpath(paths.plots_dir, "bm_state_control.pdf"), problem, bm_result)
    save_state_control_plot!(
        joinpath(paths.plots_dir, "mppi_state_control.pdf"),
        problem,
        mppi_result,
    )
    save_certification_plot!(
        joinpath(paths.plots_dir, "bm_certification.pdf"),
        problem,
        bm_result,
        bm_cert,
    )
    save_certification_plot!(
        joinpath(paths.plots_dir, "mppi_certification.pdf"),
        problem,
        mppi_result,
        mppi_cert,
    )

    write_rows_csv(joinpath(paths.data_dir, "results.csv"), rows)
    write_rows_markdown(joinpath(paths.data_dir, "results.md"), rows)

    return (
        horizon = N,
        problem = problem,
        paths = paths,
        bm = (planner = bm_result, eval = bm_eval, cert = bm_cert),
        mppi = (planner = mppi_result, eval = mppi_eval, cert = mppi_cert),
        rows = rows,
    )
end

function run_all_horizons(cfg::NonlinearBMvsMPPIConfig)
    global_paths = benchmark_output_paths(cfg)
    all_runs = [run_single_horizon(cfg, N) for N in cfg.horizons]
    summary_rows = reduce(vcat, [run.rows for run in all_runs])

    write_rows_csv(joinpath(global_paths.summary_dir, "summary.csv"), summary_rows)
    write_rows_markdown(joinpath(global_paths.summary_dir, "summary.md"), summary_rows)
    save_summary_plots!(global_paths, summary_rows)
    save_notes!(global_paths, cfg)

    return (; paths = global_paths, runs = all_runs, summary_rows)
end
