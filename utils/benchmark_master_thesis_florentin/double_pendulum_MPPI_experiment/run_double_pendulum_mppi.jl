include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import JLD2
import LinearAlgebra as LA
import MathematicalSystems as MS
import Random

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "pendulum",
        "double_pendulum.jl",
    ),
)
const DP = DoublePendulum

Base.@kwdef struct DoublePendulumMPPIConfig
    l1::Float64 = 1.0
    l2::Float64 = 1.0
    m1::Float64 = 1.0
    m2::Float64 = 1.0
    g::Float64 = 9.81
    objective::String = "benchmark_up_convex"

    Δt::Float64 = 0.05
    hx::SVector{4, Float64} = SVector(5 * (pi / 180), 5 * (pi / 180), 0.25, 0.25)
    periodic_dims::SVector{2, Int} = SVector(1, 2)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 51
    input_values::Tuple{Vararg{Float64}} = vals_tuple = Tuple(-6.5:0.25:6.5)

    terminal_radius::Float64 = 0.3
    use_terminal_john_ellipsoid::Bool = true
    terminal_john_shrink::Float64 = 1.0
    terminal_success_distance2::Float64 = 1.0
    terminal_outside_weight::Float64 = 1.0e5
    terminal_center_weight::Float64 = 50.0
    terminal_endpoint_center_weight::Float64 = 80.0
    mppi_success_criterion::Symbol = :target_set # :target_set or :terminal_ellipsoid
    mppi_input_norm_weight::Float64 = 0.08
    mppi_input_delta_weight::Float64 = 0.04
    λ::Float64 = 0.0000005
    maxδx::Float64 = 180.0
    maxδu::Float64 = 160.0
    state_scaling_mode::Symbol = :matrix # :none, :std, or :matrix
    state_scaling_std_floor::Float64 = 1.0e-3
    state_scaling_matrix::Matrix{Float64} =
        Matrix{Float64}(LA.Diagonal([1.5, 1.5, 1.5, 1.5]))
    # state_scaling_matrix::Matrix{Float64} = Matrix{Float64}(
    #     LA.Diagonal([1/pi, 1/pi, 1/6.0, 1/6.0]),
    # )
    input_scaling_mode::Symbol = :matrix # :none or :matrix
    input_scaling_matrix::Matrix{Float64} = Matrix{Float64}(LA.Diagonal([0.0001]))
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
    )
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.01, 0.01), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)
    adaptive_linearization_boxes::Bool = false
    # These radii are physical coordinates. The symbolic certifier applies
    # state_scaling internally when solving the LMI, as in the vehicle benchmarks.
    ΔX_initial::Vector{Float64} = [0.02, 0.02, 0.02, 0.02]
    ΔX_min::Vector{Float64} = [1.0e-6, 1.0e-6, 1.0e-5, 1.0e-5]
    ΔX_max::Vector{Float64} = [0.75, 0.75, 3.0, 3.0]
    ΔU_initial::Vector{Float64} = [0.02]
    ΔU_min::Vector{Float64} = [1.0e-4]
    ΔU_max::Vector{Float64} = [3.0]
    adaptive_box_growth::Float64 = 1.1
    adaptive_box_safety::Float64 = 1.01
    adaptive_box_max_iters::Int = 1
    adaptive_box_atol::Float64 = 1.0e-8
    adaptive_box_verbose::Bool = false
    adaptive_box_search_scales::Vector{Float64} = [1.0]
    adaptive_box_objective::Symbol = :first_consistent
    adaptive_box_keep_first_consistent::Bool = false

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    plot_mp4::Bool = true
    verbose::Bool = false

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 1

    planning_input_scale::Float64 = 0.75
    mppi_nsamples::Int = 1800
    mppi_niter::Int = 20
    mppi_λ::Float64 = 1.75
    mppi_noise_u::Float64 = 0.8
    rng_seed::Int = 3

    save_generated_trajectory::Bool = true
    load_saved_trajectory::Bool = true
    saved_trajectory_path::String =
        joinpath(@__DIR__, "outputs", "saved_mppi_trajectory.jld2")
    terminal_refinement_steps::Int = 2
    terminal_refinement_nsamples::Int = 6000
    terminal_refinement_niter::Int = 35
    terminal_refinement_λ::Float64 = 0.35
    terminal_refinement_noise_u::Float64 = 3.0
    terminal_refinement_center_weight::Float64 = 1500.0
    terminal_refinement_input_weight::Float64 = 0.05
    terminal_refinement_input_delta_weight::Float64 = 0.02
    terminal_refinement_rng_seed::Int = 3003
end

function build_terminal_john_ellipsoid(lower, upper; shrink = 1.0)
    lb = Float64.(collect(lower))
    ub = Float64.(collect(upper))

    length(lb) == length(ub) ||
        throw(ArgumentError("target lower/upper bounds must have the same length"))
    all(isfinite, lb) && all(isfinite, ub) ||
        throw(ArgumentError("target bounds must be finite"))
    all(ub .> lb) || throw(
        ArgumentError("each target upper bound must be strictly larger than lower bound"),
    )
    isfinite(shrink) && 0.0 < shrink <= 1.0 ||
        throw(ArgumentError("terminal_john_shrink must be finite and in (0, 1]"))

    radii = 0.5 .* (ub .- lb)
    center = 0.5 .* (lb .+ ub)
    shape = Matrix(LA.Diagonal(1.0 ./ (shrink .* radii) .^ 2))
    LA.isposdef(LA.Symmetric(shape)) ||
        throw(ArgumentError("terminal John ellipsoid shape must be positive definite"))

    return center, shape
end

function build_terminal_john_ellipsoid(target_set; shrink = 1.0)
    hasproperty(target_set, :lb) && hasproperty(target_set, :ub) ||
        throw(ArgumentError("target_set must expose lb and ub bounds"))
    return build_terminal_john_ellipsoid(target_set.lb, target_set.ub; shrink = shrink)
end

function terminal_john_ellipsoid_data(problem, cfg::DoublePendulumMPPIConfig)
    cfg.use_terminal_john_ellipsoid || return nothing
    terminal_center, terminal_shape =
        build_terminal_john_ellipsoid(problem.target_set; shrink = cfg.terminal_john_shrink)
    return (; terminal_center, terminal_shape)
end

function terminal_ellipsoidal_distance2(
    x,
    terminal_center,
    terminal_shape,
    cfg::DoublePendulumMPPIConfig,
)
    e = collect(Float64, periodic_state_error(x, terminal_center, cfg))
    return Float64(e' * terminal_shape * e)
end

function truncate_candidate_at_first_terminal_ellipsoid_hit(
    cand,
    wrap_state,
    terminal_data,
    cfg::DoublePendulumMPPIConfig;
    threshold::Float64 = 1.0,
)
    cand === nothing && return nothing
    xs = collect(ST.enum_elems(cand.x_traj))
    hit_idx = findfirst(xs) do x
        dT2 = terminal_ellipsoidal_distance2(
            wrap_state(x),
            terminal_data.terminal_center,
            terminal_data.terminal_shape,
            cfg,
        )
        return dT2 <= threshold + 1.0e-8
    end

    if hit_idx === nothing || hit_idx <= 1 || hit_idx >= length(xs)
        return cand
    end

    us = collect(ST.enum_elems(cand.u_traj))
    return OP.CandidateTrajectory(
        ST.Trajectory(xs[1:hit_idx]),
        ST.Trajectory(us[1:(hit_idx - 1)]);
        Ts = cand.Ts,
        source = cand.source,
        metadata = cand.metadata,
    )
end

function input_bounds(problem, cfg::DoublePendulumMPPIConfig)
    if hasproperty(problem.system, :U) && problem.system.U isa UT.HyperRectangle
        return Float64(problem.system.U.lb[1]), Float64(problem.system.U.ub[1])
    end
    vals = Float64.(collect(cfg.input_values))
    return minimum(vals), maximum(vals)
end

function planning_input_bounds(problem, cfg::DoublePendulumMPPIConfig)
    0.0 < cfg.planning_input_scale <= 1.0 ||
        throw(ArgumentError("planning_input_scale must be in (0, 1]"))
    umin, umax = input_bounds(problem, cfg)
    center = 0.5 * (umin + umax)
    half_width = 0.5 * (umax - umin)
    return center - cfg.planning_input_scale * half_width,
    center + cfg.planning_input_scale * half_width
end

function planning_input_domain(problem, cfg::DoublePendulumMPPIConfig)
    umin_plan, umax_plan = planning_input_bounds(problem, cfg)
    return UT.HyperRectangle([umin_plan], [umax_plan])
end

function matrix_state_scaling(M)
    S = Matrix{Float64}(M)
    size(S, 1) == size(S, 2) || error("state_scaling_matrix must be square")
    all(isfinite, S) || error("state_scaling_matrix entries must be finite")
    LA.rank(S) == size(S, 1) || error("state_scaling_matrix must be invertible")
    return S
end

function trajectory_state_std_scaling(candidate; floor::Float64 = 1.0e-3)
    xs = collect(ST.enum_elems(candidate.x_traj))
    isempty(xs) && return Matrix{Float64}(LA.I, 4, 4)

    nx = length(first(xs))
    n = length(xs)
    means = [sum(Float64(x[i]) for x in xs) / n for i in 1:nx]

    if n == 1
        return Matrix{Float64}(LA.I, nx, nx)
    end

    stds = [sqrt(sum((Float64(x[i]) - means[i])^2 for x in xs) / (n - 1)) for i in 1:nx]
    return Matrix{Float64}(LA.Diagonal(max.(stds, floor)))
end

function resolve_state_scaling(cfg::DoublePendulumMPPIConfig, candidate)
    cfg.state_scaling_mode == :none && return nothing
    cfg.state_scaling_mode == :std &&
        return trajectory_state_std_scaling(candidate; floor = cfg.state_scaling_std_floor)
    cfg.state_scaling_mode == :matrix &&
        return matrix_state_scaling(cfg.state_scaling_matrix)
    return error(
        "Unknown state_scaling_mode $(cfg.state_scaling_mode). Use :none, :std, or :matrix.",
    )
end

function matrix_input_scaling(M)
    S = Matrix{Float64}(M)
    size(S, 1) == size(S, 2) || error("input_scaling_matrix must be square")
    all(isfinite, S) || error("input_scaling_matrix entries must be finite")
    LA.rank(S) == size(S, 1) || error("input_scaling_matrix must be invertible")
    return S
end

function resolve_input_scaling(cfg::DoublePendulumMPPIConfig)
    cfg.input_scaling_mode == :none && return nothing
    cfg.input_scaling_mode == :matrix &&
        return matrix_input_scaling(cfg.input_scaling_matrix)
    return error(
        "Unknown input_scaling_mode $(cfg.input_scaling_mode). Use :none or :matrix.",
    )
end

function scaled_input_norm2(u, input_scaling)
    uv = collect(Float64, u)
    input_scaling === nothing && return sum(abs2, uv)
    z = input_scaling \ uv
    return sum(abs2, z)
end

function scaled_input_delta_norm2(u, uprev, input_scaling)
    du = collect(Float64, u) .- collect(Float64, uprev)
    input_scaling === nothing && return sum(abs2, du)
    z = input_scaling \ du
    return sum(abs2, z)
end

function candidate_reaches_target_set(cand, prob, wrap_state)
    return any(x -> (wrap_state(x) ∈ prob.target_set), ST.enum_elems(cand.x_traj))
end

function candidate_reaches_terminal_ellipsoid(
    cand,
    wrap_state,
    terminal_data,
    cfg::DoublePendulumMPPIConfig,
)
    terminal_data === nothing && return false
    return any(ST.enum_elems(cand.x_traj)) do x
        dT2 = terminal_ellipsoidal_distance2(
            wrap_state(x),
            terminal_data.terminal_center,
            terminal_data.terminal_shape,
            cfg,
        )
        return dT2 <= cfg.terminal_success_distance2 + 1.0e-8
    end
end

function mppi_candidate_success(
    cand,
    prob,
    wrap_state,
    terminal_data,
    cfg::DoublePendulumMPPIConfig,
)
    cfg.mppi_success_criterion == :target_set &&
        return candidate_reaches_target_set(cand, prob, wrap_state)
    cfg.mppi_success_criterion == :terminal_ellipsoid &&
        return candidate_reaches_terminal_ellipsoid(cand, wrap_state, terminal_data, cfg)
    return error(
        "Unknown mppi_success_criterion $(cfg.mppi_success_criterion). Use :target_set or :terminal_ellipsoid.",
    )
end

function trajectory_matrix(traj)
    elems = collect(ST.enum_elems(traj))
    isempty(elems) && return zeros(Float64, 0, 0)
    return reduce(hcat, [Float64.(collect(e)) for e in elems])
end

function svector_trajectory_from_matrix(mat::AbstractMatrix{<:Real})
    n, k = size(mat)
    return ST.Trajectory([SVector{n, Float64}(Float64.(mat[:, i])) for i in 1:k])
end

function double_pendulum_config_metadata(cfg::DoublePendulumMPPIConfig)
    return (;
        Δt = cfg.Δt,
        nstep = cfg.nstep,
        objective = cfg.objective,
        terminal_john_shrink = cfg.terminal_john_shrink,
        terminal_success_distance2 = cfg.terminal_success_distance2,
        terminal_center_weight = cfg.terminal_center_weight,
        terminal_endpoint_center_weight = cfg.terminal_endpoint_center_weight,
        mppi_input_norm_weight = cfg.mppi_input_norm_weight,
        mppi_input_delta_weight = cfg.mppi_input_delta_weight,
        planning_input_scale = cfg.planning_input_scale,
        mppi_nsamples = cfg.mppi_nsamples,
        mppi_niter = cfg.mppi_niter,
        mppi_λ = cfg.mppi_λ,
        mppi_noise_u = cfg.mppi_noise_u,
        rng_seed = cfg.rng_seed,
    )
end

function save_candidate_trajectory!(
    path::AbstractString,
    candidate,
    cfg::DoublePendulumMPPIConfig;
    generator_success::Bool,
    generator_solve_time::Float64,
    generator_diagnostics = (;),
)
    candidate === nothing && return nothing
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        x_traj = trajectory_matrix(candidate.x_traj),
        u_traj = trajectory_matrix(candidate.u_traj),
        Ts = Float64(candidate.Ts),
        source = candidate.source,
        candidate_metadata = candidate.metadata,
        generator_success = generator_success,
        generator_solve_time = generator_solve_time,
        generator_diagnostics = generator_diagnostics,
        config_metadata = double_pendulum_config_metadata(cfg),
    )
    return path
end

function load_candidate_trajectory(path::AbstractString)
    isfile(path) || error(
        "Saved trajectory not found at $(path). Set load_saved_trajectory=false " *
        "to generate a new MPPI trajectory first, or update saved_trajectory_path.",
    )
    data = JLD2.load(path)
    x_traj = svector_trajectory_from_matrix(data["x_traj"])
    u_traj = svector_trajectory_from_matrix(data["u_traj"])
    source = haskey(data, "source") ? data["source"] : :saved_mppi
    source_symbol = source isa Symbol ? source : Symbol(source)
    metadata = haskey(data, "candidate_metadata") ? data["candidate_metadata"] : (;)
    return OP.CandidateTrajectory(
        x_traj,
        u_traj;
        Ts = Float64(data["Ts"]),
        source = :saved_mppi,
        metadata = (;
            loaded_from = path,
            saved_source = source_symbol,
            saved_metadata = metadata,
        ),
    )
end

mutable struct CachedDoublePendulumTrajectoryGenerator <: OP.AbstractHeuristicGenerator
    inner::Any
    cfg::DoublePendulumMPPIConfig
    problem::Any
    candidate::Any
    seed::Any
    success::Bool
    solve_time_sec::Float64
    diagnostics::Any
end

function cached_double_pendulum_trajectory_generator(inner, cfg::DoublePendulumMPPIConfig)
    return CachedDoublePendulumTrajectoryGenerator(
        inner,
        cfg,
        nothing,
        nothing,
        nothing,
        false,
        0.0,
        (;),
    )
end

function OP.set_problem!(gen::CachedDoublePendulumTrajectoryGenerator, problem)
    gen.problem = problem
    OP.set_problem!(gen.inner, problem)
    gen.candidate = nothing
    gen.seed = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    gen.diagnostics = (;)
    return gen
end

function OP.generate!(gen::CachedDoublePendulumTrajectoryGenerator)
    gen.problem === nothing && error("Call set_problem!(gen, problem) first.")
    cfg = gen.cfg

    if cfg.load_saved_trajectory
        t = @elapsed begin
            gen.candidate = load_candidate_trajectory(cfg.saved_trajectory_path)
            gen.seed = nothing
            gen.success = gen.inner.config.success_fun(gen.problem, gen.candidate)
        end
        gen.solve_time_sec = t
        gen.diagnostics = (;
            loaded_saved_trajectory = true,
            saved_trajectory_path = cfg.saved_trajectory_path,
        )
        return gen
    end

    OP.generate!(gen.inner)
    gen.candidate = OP.get_trajectory(gen.inner)
    gen.seed = OP.get_seed(gen.inner)
    gen.success = OP.get_success(gen.inner)
    gen.solve_time_sec = OP.get_solve_time(gen.inner)

    saved_path = nothing
    if cfg.save_generated_trajectory && gen.candidate !== nothing
        saved_path = save_candidate_trajectory!(
            cfg.saved_trajectory_path,
            gen.candidate,
            cfg;
            generator_success = gen.success,
            generator_solve_time = gen.solve_time_sec,
            generator_diagnostics = gen.inner.diagnostics,
        )
    end

    gen.diagnostics = (;
        gen.inner.diagnostics...,
        loaded_saved_trajectory = false,
        saved_trajectory_path = saved_path,
    )
    return gen
end

OP.get_trajectory(gen::CachedDoublePendulumTrajectoryGenerator) = gen.candidate
OP.get_seed(gen::CachedDoublePendulumTrajectoryGenerator) = gen.seed
OP.get_success(gen::CachedDoublePendulumTrajectoryGenerator) = gen.success
OP.get_solve_time(gen::CachedDoublePendulumTrajectoryGenerator) = gen.solve_time_sec

function terminal_refinement_rollout(f_disc, wrap_state, problem, x0, u_seq)
    x = x0
    xs = Any[]
    for u in u_seq
        x = wrap_state(f_disc(x, u))
        !(x ∈ problem.system.X) && return nothing
        push!(xs, x)
    end
    return xs
end

function terminal_refinement_sequence_cost(
    f_disc,
    wrap_state,
    problem,
    x0,
    u_seq,
    previous_u,
    terminal_data,
    cfg::DoublePendulumMPPIConfig,
    input_scaling,
)
    xs = terminal_refinement_rollout(f_disc, wrap_state, problem, x0, u_seq)
    xs === nothing && return 1.0e12

    final_d2 = terminal_ellipsoidal_distance2(
        xs[end],
        terminal_data.terminal_center,
        terminal_data.terminal_shape,
        cfg,
    )
    J = cfg.terminal_refinement_center_weight * final_d2

    for k in eachindex(u_seq)
        J +=
            cfg.terminal_refinement_input_weight *
            scaled_input_norm2(u_seq[k], input_scaling)
        if k == 1
            previous_u !== nothing && (
                J +=
                    cfg.terminal_refinement_input_delta_weight *
                    scaled_input_delta_norm2(u_seq[k], previous_u, input_scaling)
            )
        else
            J +=
                cfg.terminal_refinement_input_delta_weight *
                scaled_input_delta_norm2(u_seq[k], u_seq[k - 1], input_scaling)
        end
    end
    return J
end

function _terminal_refinement_metadata(candidate, info)
    candidate.metadata isa NamedTuple &&
        return merge(candidate.metadata, (; terminal_refinement = info))
    return (; previous_metadata = candidate.metadata, terminal_refinement = info)
end

function refine_candidate_to_terminal_center(
    candidate,
    problem,
    system_cfg,
    cfg::DoublePendulumMPPIConfig,
)
    steps = max(0, cfg.terminal_refinement_steps)
    steps == 0 && return candidate, (; appended_steps = 0, reason = :disabled)

    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
    terminal_data === nothing &&
        return candidate, (; appended_steps = 0, reason = :no_terminal_ellipsoid)

    xs_existing = collect(ST.enum_elems(candidate.x_traj))
    us_existing = collect(ST.enum_elems(candidate.u_traj))
    isempty(xs_existing) &&
        return candidate, (; appended_steps = 0, reason = :empty_candidate)

    disc_system = ST.discretize_continuous_system(
        system_cfg.concrete_system,
        cfg.Δt;
        num_substeps = cfg.seed_num_substeps,
    )
    f_disc = MS.mapping(disc_system)
    wrap_state = build_pendulum_wrap_state(cfg)
    input_scaling = resolve_input_scaling(cfg)
    plan_U = planning_input_domain(problem, cfg)
    project_input = u -> project_input_to_domain(u, plan_U)

    x0 = wrap_state(xs_existing[end])
    previous_u = isempty(us_existing) ? nothing : us_existing[end]
    original_d2 = terminal_ellipsoidal_distance2(
        x0,
        terminal_data.terminal_center,
        terminal_data.terminal_shape,
        cfg,
    )

    u_nom = [project_input([0.0]) for _ in 1:steps]
    best_seq = deepcopy(u_nom)
    best_cost = terminal_refinement_sequence_cost(
        f_disc,
        wrap_state,
        problem,
        x0,
        best_seq,
        previous_u,
        terminal_data,
        cfg,
        input_scaling,
    )

    rng = Random.MersenneTwister(cfg.terminal_refinement_rng_seed)
    λ = max(cfg.terminal_refinement_λ, eps(Float64))
    nsamples = max(1, cfg.terminal_refinement_nsamples)

    for _ in 1:max(1, cfg.terminal_refinement_niter)
        eps_samples = Vector{Vector{Vector{Float64}}}(undef, nsamples)
        costs = Vector{Float64}(undef, nsamples)

        for s in 1:nsamples
            eps_seq = Vector{Vector{Float64}}(undef, steps)
            u_roll = Vector{Vector{Float64}}(undef, steps)
            for k in 1:steps
                ϵ = [cfg.terminal_refinement_noise_u * Random.randn(rng)]
                eps_seq[k] = ϵ
                u_roll[k] = project_input(u_nom[k] .+ ϵ)
            end
            eps_samples[s] = eps_seq
            costs[s] = terminal_refinement_sequence_cost(
                f_disc,
                wrap_state,
                problem,
                x0,
                u_roll,
                previous_u,
                terminal_data,
                cfg,
                input_scaling,
            )
        end

        β = minimum(costs)
        weights = exp.(-(costs .- β) ./ λ)
        weight_sum = sum(weights)
        if !isfinite(weight_sum) || weight_sum <= eps(Float64)
            weights .= 0.0
            weights[argmin(costs)] = 1.0
        else
            weights ./= weight_sum
        end

        u_new = Vector{Vector{Float64}}(undef, steps)
        for k in 1:steps
            δu = zeros(length(u_nom[k]))
            for s in 1:nsamples
                δu .+= weights[s] .* eps_samples[s][k]
            end
            u_new[k] = project_input(u_nom[k] .+ δu)
        end

        cost_new = terminal_refinement_sequence_cost(
            f_disc,
            wrap_state,
            problem,
            x0,
            u_new,
            previous_u,
            terminal_data,
            cfg,
            input_scaling,
        )
        u_nom = u_new
        if cost_new < best_cost
            best_cost = cost_new
            best_seq = deepcopy(u_new)
        end
    end

    best_prefix = 0
    best_prefix_d2 = original_d2
    best_prefix_xs = Any[]
    for h in 1:steps
        xs_h = terminal_refinement_rollout(f_disc, wrap_state, problem, x0, best_seq[1:h])
        xs_h === nothing && continue
        d2_h = terminal_ellipsoidal_distance2(
            xs_h[end],
            terminal_data.terminal_center,
            terminal_data.terminal_shape,
            cfg,
        )
        if d2_h < best_prefix_d2
            best_prefix = h
            best_prefix_d2 = d2_h
            best_prefix_xs = xs_h
        end
    end

    info = (;
        appended_steps = best_prefix,
        original_terminal_d2 = original_d2,
        refined_terminal_d2 = best_prefix_d2,
        candidate_full_horizon_cost = best_cost,
        controls = best_prefix == 0 ? Float64[] : [best_seq[k][1] for k in 1:best_prefix],
    )

    best_prefix == 0 && return candidate, info

    nu = length(best_seq[1])
    appended_us = [SVector{nu, Float64}(best_seq[k]) for k in 1:best_prefix]
    refined_candidate = OP.CandidateTrajectory(
        ST.Trajectory([xs_existing; best_prefix_xs]),
        ST.Trajectory([us_existing; appended_us]);
        Ts = candidate.Ts,
        source = :saved_mppi_terminal_refinement,
        metadata = _terminal_refinement_metadata(candidate, info),
    )
    return refined_candidate, info
end

function refine_saved_double_pendulum_trajectory_to_terminal_center!(
    cfg::DoublePendulumMPPIConfig = DoublePendulumMPPIConfig();
    save_plots::Bool = true,
)
    system_cfg = build_double_pendulum_system_cfg(cfg; pendulum_module = DP)
    control_cfg = build_double_pendulum_control_cfg(cfg; pendulum_module = DP)
    problem = control_cfg.problem
    candidate = load_candidate_trajectory(cfg.saved_trajectory_path)
    refined_candidate, info =
        refine_candidate_to_terminal_center(candidate, problem, system_cfg, cfg)

    wrap_state = build_pendulum_wrap_state(cfg)
    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
    success =
        mppi_candidate_success(refined_candidate, problem, wrap_state, terminal_data, cfg)
    save_candidate_trajectory!(
        cfg.saved_trajectory_path,
        refined_candidate,
        cfg;
        generator_success = success,
        generator_solve_time = 0.0,
        generator_diagnostics = (; terminal_refinement = info),
    )

    plot_paths = if save_plots
        save_double_pendulum_plots!(
            joinpath(cfg.output_root, cfg.plot_subdir),
            problem,
            refined_candidate;
            basename = "mppi_candidate_traj",
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            angles_title = "double_pendulum_mppi - candidate angles",
            velocities_title = "double_pendulum_mppi - candidate velocities",
            state_title = "double_pendulum_mppi - candidate states",
            control_title = "double_pendulum_mppi - candidate control",
            phase_title = "double_pendulum_mppi - candidate phase portraits",
        )
    else
        (;)
    end

    println("terminal refinement appended_steps = ", info.appended_steps)
    println(
        "terminal refinement d2: ",
        info.original_terminal_d2,
        " -> ",
        info.refined_terminal_d2,
    )
    println("saved refined trajectory: ", cfg.saved_trajectory_path)
    hasproperty(plot_paths, :velocities_path) &&
        println("saved velocity plot: ", plot_paths.velocities_path)
    return (;
        candidate = refined_candidate,
        info,
        trajectory_path = cfg.saved_trajectory_path,
        plot_paths,
    )
end

function build_double_pendulum_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::DoublePendulumMPPIConfig;
    input_mapping,
)
    seed_gen = build_centered_generator(
        problem,
        cfg;
        input_mapping = input_mapping,
        jacobian_bound = DP.jacobian_bound(;
            l1 = cfg.l1,
            l2 = cfg.l2,
            m1 = cfg.m1,
            m2 = cfg.m2,
            g = cfg.g,
            ωmax = 4.0,
        ),
        x0_provider = _ -> control_cfg.x0,
        num_substeps = cfg.seed_num_substeps,
        trajectory_mode = cfg.seed_trajectory_mode,
    )

    disc_system = ST.discretize_continuous_system(
        system_cfg.concrete_system,
        cfg.Δt;
        num_substeps = cfg.seed_num_substeps,
    )
    f_disc = MS.mapping(disc_system)
    wrap_state = build_pendulum_wrap_state(cfg)
    target_center = UT.get_center(problem.target_set)
    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
    input_scaling = resolve_input_scaling(cfg)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)
    noise_sampler = (rng, u, k) -> [cfg.mppi_noise_u * Random.randn(rng)]
    plan_U = planning_input_domain(problem, cfg)
    project_input = u -> project_input_to_domain(u, plan_U)

    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return collect(ST.enum_elems(seed_cand.x_traj))
    end

    trajectory_cost = function (prob, cand)
        xs = collect(ST.enum_elems(cand.x_traj))
        us = collect(ST.enum_elems(cand.u_traj))
        reference_states = get_reference_states()

        BAD_COST = 1.0e12
        isempty(xs) && return BAD_COST

        J = 0.0
        hit_target = false
        hit_index = length(xs)

        for k in eachindex(xs)
            xw = wrap_state(xs[k])
            !(xw ∈ prob.system.X) && return BAD_COST

            if terminal_data === nothing && xw ∈ prob.target_set
                hit_target = true
                hit_index = k
                break
            end

            e_goal = periodic_state_error(xw, target_center, cfg)
            J += 0.2
            J += 1.0 * e_goal[1]^2
            J += 3.5 * e_goal[2]^2
            J += 0.08 * e_goal[3]^2
            J += 0.08 * e_goal[4]^2

            if reference_states !== nothing
                xref = wrap_state(reference_states[min(k, length(reference_states))])
                e_ref = periodic_state_error(xw, xref, cfg)
                J += 1.2 * e_ref[1]^2
                J += 1.8 * e_ref[2]^2
                J += 0.03 * e_ref[3]^2
                J += 0.03 * e_ref[4]^2
            end
        end

        if terminal_data !== nothing
            terminal_distances = [
                terminal_ellipsoidal_distance2(
                    wrap_state(x),
                    terminal_data.terminal_center,
                    terminal_data.terminal_shape,
                    cfg,
                ) for x in xs
            ]
            hit_idx = findfirst(
                d -> d <= cfg.terminal_success_distance2 + 1.0e-8,
                terminal_distances,
            )
            hit_index = hit_idx === nothing ? argmin(terminal_distances) : hit_idx
            hit_target = hit_idx !== nothing
        end

        xT = wrap_state(xs[min(hit_index, length(xs))])
        xN = wrap_state(xs[end])
        eT = periodic_state_error(xT, target_center, cfg)
        J += 45.0 * eT[1]^2
        J += 90.0 * eT[2]^2
        J += 4.0 * eT[3]^2
        J += 4.0 * eT[4]^2

        if terminal_data !== nothing
            dT2 = terminal_ellipsoidal_distance2(
                xT,
                terminal_data.terminal_center,
                terminal_data.terminal_shape,
                cfg,
            )
            J += cfg.terminal_outside_weight * max(0.0, dT2 - 1.0)^2
            J += cfg.terminal_center_weight * dT2

            dN2 = terminal_ellipsoidal_distance2(
                xN,
                terminal_data.terminal_center,
                terminal_data.terminal_shape,
                cfg,
            )
            J += cfg.terminal_endpoint_center_weight * dN2
        end

        if !hit_target
            J += 250.0
        end

        last_u_index = min(length(us), max(hit_index - 1, 0))
        for k in 1:last_u_index
            uk = collect(Float64, us[k])
            J += cfg.mppi_input_norm_weight * scaled_input_norm2(uk, input_scaling)

            if k >= 2
                J +=
                    cfg.mppi_input_delta_weight *
                    scaled_input_delta_norm2(uk, us[k - 1], input_scaling)
            end
        end

        return J
    end

    success_fun =
        (prob, cand) -> mppi_candidate_success(cand, prob, wrap_state, terminal_data, cfg)

    postprocess_candidate = if terminal_data === nothing
        OP._default_mppi_postprocess
    else
        (_prob, _mppi_cfg, cand) -> truncate_candidate_at_first_terminal_ellipsoid_hit(
            cand,
            wrap_state,
            terminal_data,
            cfg;
            threshold = cfg.terminal_success_distance2,
        )
    end

    mppi_cfg = OP.MPPIConfig(
        cfg.Δt,
        cfg.nstep,
        cfg.mppi_nsamples,
        cfg.mppi_niter,
        cfg.mppi_λ,
        _ -> control_cfg.x0,
        discrete_dynamics,
        noise_sampler,
        project_input,
        trajectory_cost,
        success_fun,
        (prob, x) -> wrap_state(x);
        truncate_at_first_target = terminal_data === nothing,
        postprocess_candidate = postprocess_candidate,
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function _build_double_pendulum_certifier(
    problem,
    _system_cfg,
    _control_cfg,
    cfg::DoublePendulumMPPIConfig,
    state_scaling,
)
    symbolic_builder = build_double_pendulum_symbolic_builder(cfg; pendulum_module = DP)
    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        terminal_center = terminal_data === nothing ? nothing :
                          terminal_data.terminal_center,
        terminal_shape = terminal_data === nothing ? nothing : terminal_data.terminal_shape,
        state_scaling = state_scaling,
        adaptive_linearization_boxes = cfg.adaptive_linearization_boxes,
        ΔX_initial = cfg.ΔX_initial,
        ΔX_min = cfg.ΔX_min,
        ΔX_max = cfg.ΔX_max,
        ΔU_initial = cfg.ΔU_initial,
        ΔU_min = cfg.ΔU_min,
        ΔU_max = cfg.ΔU_max,
        adaptive_box_growth = cfg.adaptive_box_growth,
        adaptive_box_safety = cfg.adaptive_box_safety,
        adaptive_box_max_iters = cfg.adaptive_box_max_iters,
        adaptive_box_atol = cfg.adaptive_box_atol,
        adaptive_box_verbose = cfg.adaptive_box_verbose,
        adaptive_box_search_scales = cfg.adaptive_box_search_scales,
        adaptive_box_objective = cfg.adaptive_box_objective,
        adaptive_box_keep_first_consistent = cfg.adaptive_box_keep_first_consistent,
    )
    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        build_backend(; verbose = cfg.verbose),
        opts,
    )
    return SC.EllipsoidalBackwardCertifier(cert_cfg, symbolic_builder)
end

mutable struct DoublePendulumScalingCertifier <: SC.AbstractSymbolicCertifier
    problem::Any
    system_cfg::Any
    control_cfg::Any
    cfg::DoublePendulumMPPIConfig
    candidate::Any
    inner::Any
    result::Any
    success::Bool
    solve_time_sec::Float64
end

function build_double_pendulum_certifier(
    problem,
    system_cfg,
    control_cfg,
    cfg::DoublePendulumMPPIConfig,
)
    return DoublePendulumScalingCertifier(
        problem,
        system_cfg,
        control_cfg,
        cfg,
        nothing,
        nothing,
        nothing,
        false,
        0.0,
    )
end

function SC.set_problem!(
    cert::DoublePendulumScalingCertifier,
    problem::Dionysos.Problem.ProblemType,
)
    cert.problem = problem
    cert.candidate = nothing
    cert.inner = nothing
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function SC.set_trajectory!(cert::DoublePendulumScalingCertifier, candidate)
    cert.candidate = candidate
    cert.inner = nothing
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function SC.certify!(cert::DoublePendulumScalingCertifier)
    state_scaling = resolve_state_scaling(cert.cfg, cert.candidate)
    cert.inner = _build_double_pendulum_certifier(
        cert.problem,
        cert.system_cfg,
        cert.control_cfg,
        cert.cfg,
        state_scaling,
    )
    SC.set_problem!(cert.inner, cert.problem)
    SC.set_trajectory!(cert.inner, cert.candidate)
    SC.certify!(cert.inner)

    cert.result = SC.get_result(cert.inner)
    cert.success = SC.get_success(cert.inner)
    cert.solve_time_sec = SC.get_solve_time(cert.inner)
    return cert
end

SC.get_result(cert::DoublePendulumScalingCertifier) = cert.result
SC.get_success(cert::DoublePendulumScalingCertifier) = cert.success
SC.get_solve_time(cert::DoublePendulumScalingCertifier) = cert.solve_time_sec

function active_state_scaling(cert::DoublePendulumScalingCertifier)
    cert.inner === nothing && return nothing
    return cert.inner.config.options.state_scaling
end

function main(cfg::DoublePendulumMPPIConfig = DoublePendulumMPPIConfig())
    Random.seed!(cfg.rng_seed)

    input_mapping = build_double_pendulum_input_mapping(cfg)

    generator_builder =
        (problem, system_cfg, control_cfg, cfg_) -> begin
            inner_gen = build_double_pendulum_mppi_generator(
                problem,
                system_cfg,
                control_cfg,
                cfg_;
                input_mapping = input_mapping,
            )
            return cached_double_pendulum_trajectory_generator(inner_gen, cfg_)
        end

    certifier_builder =
        (problem, system_cfg, control_cfg, cfg_) ->
            build_double_pendulum_certifier(problem, system_cfg, control_cfg, cfg_)

    save_artifacts! = function (run_result)
        save_double_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = "ellipsoids",
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            angles_title = "double_pendulum_mppi - certification angles",
            velocities_title = "double_pendulum_mppi - certification velocities",
            state_title = "double_pendulum_mppi - certification states",
            control_title = "double_pendulum_mppi - certification control",
            phase_title = "double_pendulum_mppi - certification phase portraits",
        )
        animation_paths = save_double_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = "ellipsoids",
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            l1 = cfg.l1,
            l2 = cfg.l2,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "double_pendulum_mppi - certification",
        )
        return (; gif = animation_paths.gif_path, mp4 = animation_paths.mp4_path)
    end

    run_result = run_benchmark(
        cfg;
        scenario_name = "double_pendulum_mppi",
        build_concrete_system = () ->
            build_double_pendulum_system_cfg(cfg; pendulum_module = DP),
        build_control_problem = () ->
            build_double_pendulum_control_cfg(cfg; pendulum_module = DP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )

    warmup_candidate = OP.get_seed(run_result.solver.generator)
    if warmup_candidate !== nothing
        save_double_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            warmup_candidate;
            basename = "abstract_traj_used_as_warmup",
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            angles_title = "double_pendulum_mppi - warmup angles",
            velocities_title = "double_pendulum_mppi - warmup velocities",
            state_title = "double_pendulum_mppi - warmup states",
            control_title = "double_pendulum_mppi - warmup control",
            phase_title = "double_pendulum_mppi - warmup phase portraits",
        )
        save_double_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            warmup_candidate;
            basename = "abstract_traj_used_as_warmup",
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            l1 = cfg.l1,
            l2 = cfg.l2,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "double_pendulum_mppi - warmup",
        )
    else
        println("warmup trajectory unavailable; skipping warmup plots.")
    end

    save_double_pendulum_plots!(
        run_result.outputs.plots_dir,
        run_result.problem,
        run_result.nominal_candidate;
        basename = "mppi_candidate_traj",
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        angles_title = "double_pendulum_mppi - candidate angles",
        velocities_title = "double_pendulum_mppi - candidate velocities",
        state_title = "double_pendulum_mppi - candidate states",
        control_title = "double_pendulum_mppi - candidate control",
        phase_title = "double_pendulum_mppi - candidate phase portraits",
    )
    save_double_pendulum_animation!(
        run_result.outputs.animations_dir,
        run_result.problem,
        run_result.nominal_candidate;
        basename = "mppi_candidate_traj",
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        l1 = cfg.l1,
        l2 = cfg.l2,
        save_gif = cfg.plot_gif,
        save_mp4 = cfg.plot_mp4,
        title = "double_pendulum_mppi - candidate",
    )

    stat_result = run_kappa_statistical_check(run_result; n_samples = 500)
    save_kappa_statistical_plots!(stat_result; wrap_angles = true)

    return run_result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
