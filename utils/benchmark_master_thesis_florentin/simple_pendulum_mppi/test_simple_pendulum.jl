include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
using LaTeXStrings
import LinearAlgebra as LA
import Printf
import MathematicalSystems as MS
import Random
import Statistics

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pendulum", "simple_pendulum.jl"))
const SP = SimplePendulum

Base.@kwdef struct SimplePendulumMPPIConfig
    l::Float64 = 1.0
    g::Float64 = 9.81
    objective::String = "benchmark_up_convex"

    Δt::Float64 = 0.1
    hx::SVector{2, Float64} = SVector(5 * (pi / 180), 0.25)
    periodic_dims::SVector{1, Int} = SVector(1)
    periodic_periods::SVector{1, Float64} = SVector(2pi)
    periodic_start::SVector{1, Float64} = SVector(-pi)
    nstep::Int = 50
    input_values::Tuple{Vararg{Float64}} = Tuple(-2.6:0.25:2.6)

    terminal_radius::Float64 = 0.25
    use_terminal_inner_ellipsoid::Bool = true
    terminal_shrink::Float64 = 0.85
    terminal_outside_weight::Float64 = 1.0e6
    terminal_center_weight::Float64 = 1.0e5
    terminal_success_distance2::Float64 = 0.1
    velocity_running_weight::Float64 = 0.15
    input_margin_weight::Float64 = 0.0
    terminal_velocity_weight::Float64 = 8.0
    λ::Float64 = 0.005
    maxδx::Float64 = 80.0
    maxδu::Float64 = 40.0
    use_state_scaling::Bool = true
    state_scaling::Vector{Float64} = [0.85 * 15.0 * pi / 180.0, 0.85]
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.0002, 0.0002),
        IA.interval(-0.0002, 0.0002),
    )
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.0002, 0.0002), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)
    adaptive_linearization_boxes::Bool = true
    ΔX_initial::Vector{Float64} = [0.02, 0.02]
    ΔX_min::Vector{Float64} = [1.0e-3, 1.0e-3]
    ΔX_max::Vector{Float64} = [0.30, 2.00]
    ΔU_initial::Vector{Float64} = [0.02]
    ΔU_min::Vector{Float64} = [1.0e-3]
    ΔU_max::Vector{Float64} = [1.00]
    adaptive_box_growth::Float64 = 1.5
    adaptive_box_safety::Float64 = 1.15
    adaptive_box_max_iters::Int = 8
    adaptive_box_atol::Float64 = 1.0e-8
    adaptive_box_verbose::Bool = true
    adaptive_box_search_scales::Vector{Float64} = [0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 2.0]
    adaptive_box_objective::Symbol = :max_volume
    adaptive_box_keep_first_consistent::Bool = false

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    plot_mp4::Bool = true
    verbose::Bool = false

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 1

    mppi_nsamples::Int = 7800
    mppi_niter::Int = 20
    mppi_λ::Float64 = 1.75
    mppi_noise_u::Float64 = 1.8
    rng_seed::Int = 1
    kappa_statistical_samples::Int = 12000
end

function build_inner_terminal_ellipsoid(lower, upper; shrink = 0.95)
    lb = Float64.(collect(lower))
    ub = Float64.(collect(upper))

    length(lb) == length(ub) ||
        throw(ArgumentError("target lower/upper bounds must have the same length"))
    all(isfinite, lb) && all(isfinite, ub) ||
        throw(ArgumentError("target bounds must be finite"))
    all(ub .> lb) ||
        throw(ArgumentError("each target upper bound must be strictly larger than lower bound"))
    isfinite(shrink) && 0.0 < shrink <= 1.0 ||
        throw(ArgumentError("terminal_shrink must be finite and in (0, 1]"))

    radii = 0.5 .* (ub .- lb)
    all(radii .> 0.0) ||
        throw(ArgumentError("target radii must be positive"))

    terminal_center = 0.5 .* (lb .+ ub)
    terminal_shape = Matrix(LA.Diagonal(1.0 ./ (shrink .* radii) .^ 2))
    LA.isposdef(LA.Symmetric(terminal_shape)) ||
        throw(ArgumentError("terminal_shape must be positive definite"))

    return terminal_center, terminal_shape
end

function build_inner_terminal_ellipsoid(target_set; shrink = 0.85)
    hasproperty(target_set, :lb) && hasproperty(target_set, :ub) ||
        throw(ArgumentError("target_set must expose lb and ub bounds"))
    return build_inner_terminal_ellipsoid(target_set.lb, target_set.ub; shrink = shrink)
end

function check_ellipsoid_inside_box(c, P, lower, upper; atol = 1.0e-8)
    center = Float64.(collect(c))
    shape = Matrix{Float64}(P)
    lb = Float64.(collect(lower))
    ub = Float64.(collect(upper))

    length(center) == length(lb) == length(ub) ||
        throw(ArgumentError("ellipsoid center and box bounds must have the same length"))
    size(shape) == (length(center), length(center)) ||
        throw(ArgumentError("ellipsoid shape has incompatible size $(size(shape))"))
    LA.isposdef(LA.Symmetric(shape)) ||
        throw(ArgumentError("ellipsoid shape must be positive definite"))

    Q = inv(LA.Symmetric(shape))
    radii = sqrt.(max.(0.0, LA.diag(Q)))
    return all(center .+ radii .<= ub .+ atol) &&
           all(center .- radii .>= lb .- atol)
end

terminal_ellipsoidal_distance2(x, terminal_center, terminal_shape) = begin
    e = Float64.(collect(x)) .- terminal_center
    return Float64(e' * terminal_shape * e)
end

function truncate_candidate_at_first_terminal_ellipsoid_hit(
    cand,
    wrap_state,
    terminal_data;
    threshold::Float64 = 1.0,
)
    cand === nothing && return nothing
    xs = collect(ST.enum_elems(cand.x_traj))
    hit_idx = findfirst(xs) do x
        dT2 = terminal_ellipsoidal_distance2(
            wrap_state(x),
            terminal_data.terminal_center,
            terminal_data.terminal_shape,
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

function terminal_inner_ellipsoid_data(problem, cfg::SimplePendulumMPPIConfig)
    cfg.use_terminal_inner_ellipsoid || return nothing
    terminal_center, terminal_shape = build_inner_terminal_ellipsoid(
        problem.target_set;
        shrink = cfg.terminal_shrink,
    )
    inside = check_ellipsoid_inside_box(
        terminal_center,
        terminal_shape,
        problem.target_set.lb,
        problem.target_set.ub,
    )
    return (; terminal_center, terminal_shape, inside)
end

function print_terminal_inner_ellipsoid_diagnostics(problem, cfg, terminal_data)
    println("=== Terminal inner ellipsoid diagnostics ===")
    println("use_terminal_inner_ellipsoid: ", cfg.use_terminal_inner_ellipsoid)
    println("terminal_shrink: ", cfg.terminal_shrink)
    if terminal_data === nothing
        println("terminal_center: nothing")
        println("terminal_shape: nothing")
        println("terminal ellipsoid inside target set: false")
        return nothing
    end
    println("terminal_center: ", terminal_data.terminal_center)
    println("terminal_shape: ", terminal_data.terminal_shape)
    println("terminal ellipsoid inside target set: ", terminal_data.inside)
    println("terminal outside weight: ", cfg.terminal_outside_weight)
    println("terminal center weight: ", cfg.terminal_center_weight)
    println("terminal success distance dT2 threshold: ", cfg.terminal_success_distance2)
    return nothing
end

function build_simple_pendulum_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::SimplePendulumMPPIConfig;
    input_mapping,
)
    seed_gen = build_centered_generator(
        problem,
        cfg;
        input_mapping = input_mapping,
        jacobian_bound = SP.jacobian_bound(cfg.l, cfg.g),
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
    terminal_data = terminal_inner_ellipsoid_data(problem, cfg)

    discrete_dynamics = (_prob, x, u, _k, _Δt) -> f_disc(x, u)
    noise_sampler = (rng, _u, _k) -> [cfg.mppi_noise_u * Random.randn(rng)]
    project_input = u -> project_pendulum_input_to_domain(u, problem.system.U)

    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return collect(ST.enum_elems(seed_cand.x_traj))
    end

    trajectory_cost = function (prob, cand)
        xs = collect(ST.enum_elems(cand.x_traj))
        us = collect(ST.enum_elems(cand.u_traj))
        reference_states = get_reference_states()

        bad_cost = 1.0e12
        isempty(xs) && return bad_cost

        J = 0.0
        hit_target = false
        hit_index = length(xs)

        for k in eachindex(xs)
            xw = wrap_state(xs[k])
            !(xw ∈ prob.system.X) && return bad_cost

            if terminal_data === nothing && xw ∈ prob.target_set
                hit_target = true
                hit_index = k
                break
            end

            e_goal = pendulum_state_error(xw, target_center, cfg)
            J += 0.12
            J += 3.0 * e_goal[1]^2
            J += cfg.velocity_running_weight * e_goal[2]^2

            if reference_states !== nothing
                xref = wrap_state(reference_states[min(k, length(reference_states))])
                e_ref = pendulum_state_error(xw, xref, cfg)
                J += 1.0 * e_ref[1]^2
                J += 0.04 * e_ref[2]^2
            end
        end

        if terminal_data !== nothing
            terminal_distances = [
                terminal_ellipsoidal_distance2(
                    wrap_state(x),
                    terminal_data.terminal_center,
                    terminal_data.terminal_shape,
                ) for x in xs
            ]
            hit_idx = findfirst(d -> d <= cfg.terminal_success_distance2 + 1.0e-8, terminal_distances)
            hit_index = hit_idx === nothing ? argmin(terminal_distances) : hit_idx
            hit_target = hit_idx !== nothing
        end

        xT = wrap_state(xs[min(hit_index, length(xs))])
        eT = pendulum_state_error(xT, target_center, cfg)
        J += 80.0 * eT[1]^2
        J += cfg.terminal_velocity_weight * eT[2]^2

        if terminal_data !== nothing
            dT2 = terminal_ellipsoidal_distance2(
                xT,
                terminal_data.terminal_center,
                terminal_data.terminal_shape,
            )
            J += cfg.terminal_outside_weight * max(0.0, dT2 - 1.0)^2
            J += cfg.terminal_center_weight * dT2
        end

        if !hit_target
            J += 180.0
        end

        last_u_index = min(length(us), max(hit_index - 1, 0))
        umax = if hasproperty(prob.system, :U) && prob.system.U isa UT.HyperRectangle
            maximum(abs.([prob.system.U.lb[1], prob.system.U.ub[1]]))
        else
            maximum(abs.(cfg.input_values))
        end
        for k in 1:last_u_index
            uk = collect(Float64, us[k])
            J += 0.02 * uk[1]^2
            if cfg.input_margin_weight > 0.0 && umax > 0.0
                J += cfg.input_margin_weight * (abs(uk[1]) / umax)^4
            end

            if k >= 2
                duk = uk[1] - us[k - 1][1]
                J += 0.04 * duk^2
            end
        end

        return J
    end

    success_fun = if terminal_data === nothing
        (prob, cand) -> any(x -> (wrap_state(x) ∈ prob.target_set), ST.enum_elems(cand.x_traj))
    else
        (_prob, cand) -> any(ST.enum_elems(cand.x_traj)) do x
            dT2 = terminal_ellipsoidal_distance2(
                wrap_state(x),
                terminal_data.terminal_center,
                terminal_data.terminal_shape,
            )
            return dT2 <= cfg.terminal_success_distance2 + 1.0e-8
        end
    end

    postprocess_candidate = if terminal_data === nothing
        OP._default_mppi_postprocess
    else
        (_prob, _mppi_cfg, cand) ->
            truncate_candidate_at_first_terminal_ellipsoid_hit(
                cand,
                wrap_state,
                terminal_data;
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
        (_prob, x) -> wrap_state(x);
        truncate_at_first_target = terminal_data === nothing,
        postprocess_candidate = postprocess_candidate,
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function build_simple_pendulum_certifier(
    problem,
    _system_cfg,
    _control_cfg,
    cfg::SimplePendulumMPPIConfig,
)
    symbolic_builder = build_pendulum_symbolic_builder(cfg; pendulum_module = SP)
    terminal_data = terminal_inner_ellipsoid_data(problem, cfg)
    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        terminal_center = terminal_data === nothing ? nothing : terminal_data.terminal_center,
        terminal_shape = terminal_data === nothing ? nothing : terminal_data.terminal_shape,
        state_scaling = certifier_state_scaling(cfg),
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

certifier_state_scaling(cfg::SimplePendulumMPPIConfig) =
    cfg.use_state_scaling ? Float64.(cfg.state_scaling) : nothing

function _step_field(step, name::Symbol, default = nothing)
    return hasproperty(step.summary, name) ? getproperty(step.summary, name) : default
end

function _ellipsoid_volume_or_missing(E)
    E === nothing && return missing
    try
        return UT.get_volume(E)
    catch
        return missing
    end
end

function _config_namedtuple(cfg::SimplePendulumMPPIConfig)
    names = fieldnames(typeof(cfg))
    return NamedTuple{names}(getfield(cfg, name) for name in names)
end

function _with_config(cfg::SimplePendulumMPPIConfig; kwargs...)
    return SimplePendulumMPPIConfig(; merge(_config_namedtuple(cfg), NamedTuple(kwargs))...)
end

function unwrap_pendulum_angles(xs, cfg::SimplePendulumMPPIConfig)
    isempty(xs) && return Float64[]
    period = Float64(cfg.periodic_periods[1])
    angles = [Float64(collect(x)[1]) for x in xs]
    unwrapped = similar(angles)
    unwrapped[1] = angles[1]
    half_period = 0.5 * period
    for i in 2:length(angles)
        δ = angles[i] - angles[i - 1]
        while δ > half_period
            δ -= period
        end
        while δ < -half_period
            δ += period
        end
        unwrapped[i] = unwrapped[i - 1] + δ
    end
    return unwrapped
end

function trajectory_std_scaling(xs, cfg::SimplePendulumMPPIConfig; min_scale = 1.0e-2)
    length(xs) >= 2 || return max.(Float64.(cfg.state_scaling), min_scale)
    θ = unwrap_pendulum_angles(xs, cfg)
    ω = [Float64(collect(x)[2]) for x in xs]
    return max.([Statistics.std(θ), Statistics.std(ω)], min_scale)
end

function compute_trajectory_std_scaling_from_candidate(
    cand,
    cfg::SimplePendulumMPPIConfig;
    min_scale = 1.0e-2,
)
    cand === nothing && throw(ArgumentError("candidate trajectory is missing"))
    xs = collect(ST.enum_elems(cand.x_traj))
    return trajectory_std_scaling(xs, cfg; min_scale = min_scale)
end

function trajectory_range_scaling(xs, cfg::SimplePendulumMPPIConfig; min_scale = 1.0e-2)
    isempty(xs) && return max.(Float64.(cfg.state_scaling), min_scale)
    θ = unwrap_pendulum_angles(xs, cfg)
    ω = [Float64(collect(x)[2]) for x in xs]
    return max.([maximum(θ) - minimum(θ), maximum(ω) - minimum(ω)], min_scale)
end

function scaling_candidates_from_nominal(run_result, cfg::SimplePendulumMPPIConfig)
    xs = collect(ST.enum_elems(run_result.certification_candidate.x_traj))
    current = [0.85 * 15.0 * pi / 180.0, 0.85]
    candidates = [
        (; name = "no_scaling", use_state_scaling = false, scaling = nothing),
        (; name = "current", use_state_scaling = true, scaling = current),
        (; name = "half_current", use_state_scaling = true, scaling = 0.5 .* current),
        (; name = "double_current", use_state_scaling = true, scaling = 2.0 .* current),
        (;
            name = "trajectory_std",
            use_state_scaling = true,
            scaling = trajectory_std_scaling(xs, cfg),
        ),
        (;
            name = "trajectory_range",
            use_state_scaling = true,
            scaling = trajectory_range_scaling(xs, cfg),
        ),
    ]
    return candidates
end

function _finite_values(values)
    out = Float64[]
    for value in values
        value === missing && continue
        value === nothing && continue
        isfinite(Float64(value)) || continue
        push!(out, Float64(value))
    end
    return out
end

function _safe_minimum(values)
    vals = _finite_values(values)
    return isempty(vals) ? missing : minimum(vals)
end

function _safe_maximum(values)
    vals = _finite_values(values)
    return isempty(vals) ? missing : maximum(vals)
end

function _safe_mean(values)
    vals = _finite_values(values)
    return isempty(vals) ? missing : Statistics.mean(vals)
end

function _safe_median(values)
    vals = _finite_values(values)
    return isempty(vals) ? missing : Statistics.median(vals)
end

function _safe_ratio_margin(numer, denom)
    numer === nothing && return missing
    denom === nothing && return missing
    n = Float64.(collect(numer))
    d = Float64.(collect(denom))
    length(n) == length(d) || return missing
    isempty(n) && return missing
    margins = Float64[]
    for (a, b) in zip(n, d)
        if abs(b) <= eps(Float64)
            push!(margins, a >= -eps(Float64) ? Inf : -Inf)
        else
            push!(margins, a / b)
        end
    end
    return minimum(margins)
end

function _certification_volumes(cert)
    cert === nothing && return Float64[]
    volumes = Float64[]
    if hasproperty(cert, :lmi_data) &&
       cert.lmi_data !== nothing &&
       hasproperty(cert.lmi_data, :ellipsoids)
        for E in cert.lmi_data.ellipsoids
            v = _ellipsoid_volume_or_missing(E)
            v === missing || push!(volumes, Float64(v))
        end
    end
    if isempty(volumes) && hasproperty(cert, :steps)
        for step in cert.steps
            step.status == :ok || continue
            v = _ellipsoid_volume_or_missing(step.ellipsoid)
            v === missing || push!(volumes, Float64(v))
        end
    end
    return volumes
end

function _certification_volume_summary(cert)
    volumes = _certification_volumes(cert)
    return (;
        initial = cert.success && !isempty(volumes) ? last(volumes) : missing,
        minimum = isempty(volumes) ? missing : minimum(volumes),
        median = isempty(volumes) ? missing : Statistics.median(volumes),
        terminal = isempty(volumes) ? missing : first(volumes),
    )
end

function _adaptive_box_margins(cert)
    x_margins = Float64[]
    u_margins = Float64[]
    statuses = Symbol[]
    iters = Float64[]
    candidates = Float64[]
    cert === nothing && return (; x_margins, u_margins, statuses, iters, candidates)
    for step in cert.steps
        step.status == :ok || continue
        push!(statuses, _step_field(step, :adaptive_box_status, :missing))
        xm = _safe_ratio_margin(
            _step_field(step, :Xbar_radius),
            _step_field(step, :required_X_radius),
        )
        um = _safe_ratio_margin(
            _step_field(step, :Ubar_radius),
            _step_field(step, :required_U_radius),
        )
        xm === missing || push!(x_margins, Float64(xm))
        um === missing || push!(u_margins, Float64(um))
        iter = _step_field(step, :adaptive_box_iters)
        nbox = _step_field(step, :number_of_candidate_boxes)
        iter === nothing || iter === missing || push!(iters, Float64(iter))
        nbox === nothing || nbox === missing || push!(candidates, Float64(nbox))
    end
    return (; x_margins, u_margins, statuses, iters, candidates)
end

function _certified_ellipsoid_chain(cert)
    cert === nothing && return NamedTuple[]
    chain = NamedTuple[]
    for step in cert.steps
        step.status == :ok || continue
        step.ellipsoid === nothing && continue
        push!(chain, (; k = step.k, E = step.ellipsoid))
    end
    if hasproperty(cert, :lmi_data) &&
       cert.lmi_data !== nothing &&
       hasproperty(cert.lmi_data, :ellipsoids) &&
       !isempty(cert.lmi_data.ellipsoids)
        k_terminal = isempty(chain) ? 1 : maximum(r.k for r in chain) + 1
        push!(chain, (; k = k_terminal, E = first(cert.lmi_data.ellipsoids)))
    end
    return sort(chain; by = r -> r.k)
end

function _ellipsoid_principal_radii(E)
    λ = LA.eigvals(LA.Symmetric(Matrix{Float64}(E.P)))
    return sort(1.0 ./ sqrt.(max.(Float64.(λ), eps(Float64))); rev = true)
end

function _step_summary_series(cert, field::Symbol)
    ks = Int[]
    values = Any[]
    cert === nothing && return (; ks, values)
    for step in sort(cert.steps; by = s -> s.k)
        step.status == :ok || continue
        push!(ks, step.k)
        push!(values, _step_field(step, field))
    end
    return (; ks, values)
end

function _min_ratio_series(cert, numerator_field::Symbol, denominator_field::Symbol)
    ks = Int[]
    vals = Float64[]
    cert === nothing && return (; ks, vals)
    for step in sort(cert.steps; by = s -> s.k)
        step.status == :ok || continue
        ratio = _safe_ratio_margin(_step_field(step, numerator_field), _step_field(step, denominator_field))
        ratio === missing && continue
        push!(ks, step.k)
        push!(vals, Float64(ratio))
    end
    return (; ks, vals)
end

function _vec_step_series(cert, field::Symbol)
    ks = Int[]
    vals = Vector{Float64}[]
    cert === nothing && return (; ks, vals)
    for step in sort(cert.steps; by = s -> s.k)
        step.status == :ok || continue
        value = _step_field(step, field)
        value === nothing && continue
        value === missing && continue
        push!(ks, step.k)
        push!(vals, Float64.(collect(value)))
    end
    return (; ks, vals)
end

function _input_bounds(problem, cfg::SimplePendulumMPPIConfig)
    if hasproperty(problem.system, :U) && problem.system.U isa UT.HyperRectangle
        return Float64(problem.system.U.lb[1]), Float64(problem.system.U.ub[1])
    end
    vals = Float64.(collect(cfg.input_values))
    return minimum(vals), maximum(vals)
end

function _nominal_state_control_series(run_result)
    xs = [Float64.(collect(x)) for x in ST.enum_elems(run_result.certification_candidate.x_traj)]
    us = [Float64.(collect(u)) for u in ST.enum_elems(run_result.certification_candidate.u_traj)]
    return (; xs, us, ks_x = collect(1:length(xs)), ks_u = collect(1:length(us)))
end

function _savefig_pdf(fig, output_dir, filename)
    mkpath(output_dir)
    path = joinpath(output_dir, filename)
    savefig(fig, path)
    return path
end

function _plot_selected_ellipsoid_window!(
    fig,
    chain,
    xs,
    selected_k::Int;
    window::Int = 4,
    title::AbstractString = "",
)
    local_chain = [r for r in chain if abs(r.k - selected_k) <= window]
    isempty(local_chain) && (local_chain = [chain[argmin(abs.(getfield.(chain, :k) .- selected_k))]])
    kmin = max(1, minimum(r.k for r in local_chain) - 1)
    kmax = min(length(xs), maximum(r.k for r in local_chain) + 1)
    plot!(
        fig,
        [xs[i][1] for i in kmin:kmax],
        [xs[i][2] for i in kmin:kmax];
        color = THESIS_TRAJ_COLOR,
        lw = 1.1,
        marker = :circle,
        ms = 2.2,
        label = "local nominal trajectory",
    )
    used = false
    for rec in local_chain
        plot!(
            fig,
            project_ellipsoid_2d(rec.E; dims = (1, 2));
            color = THESIS_ELLIPSOID_COLOR,
            opacity = rec.k == selected_k ? 0.22 : 0.10,
            lw = rec.k == selected_k ? 1.3 : 0.7,
            label = used ? "" : "certified ellipsoids",
        )
        used = true
    end
    return fig
end

function _zoom_limits_for_window(chain, xs, selected_k; window::Int = 4)
    local_chain = [r for r in chain if abs(r.k - selected_k) <= window]
    isempty(local_chain) && (local_chain = chain)
    kmin = max(1, minimum(r.k for r in local_chain) - 1)
    kmax = min(length(xs), maximum(r.k for r in local_chain) + 1)
    bounds = Any[_trajectory_bounds_2d(xs[kmin:kmax], (1, 2))]
    append!(bounds, [_projected_ellipsoid_bounds(r.E, (1, 2)) for r in local_chain])
    return _expand_bounds_2d(_merge_bounds_2d(bounds); margin_ratio = 0.18, min_width = 0.08, min_height = 0.08)
end

function _build_zoom_phase_plot(run_result, selected_k::Int; window::Int = 4, title::AbstractString = "")
    chain = _certified_ellipsoid_chain(run_result.result.certification)
    xs = [Float64.(collect(x)) for x in ST.enum_elems(run_result.certification_candidate.x_traj)]
    lims = _zoom_limits_for_window(chain, xs, selected_k; window = window)
    fig = plot(;
        thesis_plot_kwargs(; legend = :topright, size = (640, 460))...,
        xlims = lims.xlims,
        ylims = lims.ylims,
        xlabel = L"\theta\,[\mathrm{rad}]",
        ylabel = L"\dot{\theta}\,[\mathrm{rad/s}]",
        title = title,
    )
    if hasproperty(run_result.problem, :initial_set)
        plot!(fig, run_result.problem.initial_set; dims = [1, 2], color = :seagreen, opacity = 0.12, lw = 0.7, label = L"\mathcal{X}_I")
    end
    if hasproperty(run_result.problem, :target_set)
        plot!(fig, run_result.problem.target_set; dims = [1, 2], color = :steelblue, opacity = 0.14, lw = 0.7, label = L"\mathcal{X}_T")
    end
    _plot_selected_ellipsoid_window!(fig, chain, xs, selected_k; window = window, title = title)
    return fig
end

function _build_selected_steps_zoom_plot(run_result, selected_ks)
    figs = [
        _build_zoom_phase_plot(run_result, k; window = 3, title = "k=$(k)") for k in selected_ks
    ]
    return plot(figs...; layout = (2, 2), size = (1100, 850))
end

function save_simple_pendulum_certification_diagnostics!(
    run_result,
    cfg::SimplePendulumMPPIConfig;
    basename::AbstractString = "trajectory_std_diagnostics",
)
    cert = run_result.result.certification
    plots_dir = run_result.outputs.plots_dir
    chain = _certified_ellipsoid_chain(cert)
    isempty(chain) && return (;)

    ks = [r.k for r in chain]
    volumes = [_ellipsoid_volume_or_missing(r.E) for r in chain]
    finite_volumes = Float64[v for v in volumes if v !== missing]
    min_idx = argmin([v === missing ? Inf : Float64(v) for v in volumes])
    min_k = ks[min_idx]
    initial_k = minimum(ks)
    terminal_k = maximum(ks)

    fig_vol = plot(
        ks,
        [v === missing ? NaN : Float64(v) for v in volumes];
        yscale = :log10,
        marker = :circle,
        xlabel = "k",
        ylabel = "ellipsoid volume",
        label = "volume",
        thesis_plot_kwargs(; legend = :topright, size = (700, 420))...,
    )
    volume_path = _savefig_pdf(fig_vol, plots_dir, "trajectory_std_volume_vs_step.pdf")

    radii = [_ellipsoid_principal_radii(r.E) for r in chain]
    fig_pr = plot(;
        xlabel = "k",
        ylabel = "principal radius",
        yscale = :log10,
        thesis_plot_kwargs(; legend = :topright, size = (700, 420))...,
    )
    plot!(fig_pr, ks, [r[1] for r in radii]; marker = :circle, label = "radius 1")
    plot!(fig_pr, ks, [r[2] for r in radii]; marker = :diamond, label = "radius 2")
    principal_radii_path = _savefig_pdf(fig_pr, plots_dir, "trajectory_std_principal_radii_vs_step.pdf")

    Ls = _vec_step_series(cert, :L)
    fig_L = plot(;
        xlabel = "k",
        ylabel = "L",
        yscale = :log10,
        thesis_plot_kwargs(; legend = :topright, size = (700, 420))...,
    )
    if !isempty(Ls.vals)
        nL = maximum(length.(Ls.vals))
        for i in 1:nL
            plot!(fig_L, Ls.ks, [i <= length(v) ? max(v[i], eps(Float64)) : NaN for v in Ls.vals]; marker = :circle, label = "L[$i]")
        end
    end
    remainder_path = _savefig_pdf(fig_L, plots_dir, "trajectory_std_remainder_L_vs_step.pdf")

    xm = _min_ratio_series(cert, :Xbar_radius, :required_X_radius)
    um = _min_ratio_series(cert, :Ubar_radius, :required_U_radius)
    fig_m = plot(;
        xlabel = "k",
        ylabel = "containment margin",
        thesis_plot_kwargs(; legend = :topright, size = (700, 420))...,
    )
    plot!(fig_m, xm.ks, xm.vals; marker = :circle, label = "min X box / required")
    plot!(fig_m, um.ks, um.vals; marker = :diamond, label = "min U box / required")
    hline!(fig_m, [1.0]; color = :black, ls = :dash, label = "soundness threshold")
    margins_path = _savefig_pdf(fig_m, plots_dir, "trajectory_std_box_margins_vs_step.pdf")

    Xbar = _vec_step_series(cert, :Xbar_radius)
    Xreq = _vec_step_series(cert, :required_X_radius)
    Ubar = _vec_step_series(cert, :Ubar_radius)
    Ureq = _vec_step_series(cert, :required_U_radius)
    fig_box = plot(;
        xlabel = "k",
        ylabel = "radius",
        yscale = :log10,
        thesis_plot_kwargs(; legend = :topright, size = (900, 600))...,
    )
    for i in 1:2
        plot!(fig_box, Xbar.ks, [v[i] for v in Xbar.vals]; marker = :circle, label = "Xbar[$i]")
        plot!(fig_box, Xreq.ks, [v[i] for v in Xreq.vals]; marker = :diamond, ls = :dash, label = "required X[$i]")
    end
    plot!(fig_box, Ubar.ks, [v[1] for v in Ubar.vals]; marker = :utriangle, label = "Ubar")
    plot!(fig_box, Ureq.ks, [v[1] for v in Ureq.vals]; marker = :dtriangle, ls = :dash, label = "required U")
    box_radii_path = _savefig_pdf(fig_box, plots_dir, "trajectory_std_box_radii_vs_step.pdf")

    series = _nominal_state_control_series(run_result)
    fig_nom = plot(;
        layout = (3, 1),
        thesis_plot_kwargs(; legend = :topright, size = (850, 850))...,
    )
    plot!(fig_nom[1], series.ks_x, [x[1] for x in series.xs]; marker = :circle, label = L"\theta", ylabel = L"\theta")
    plot!(fig_nom[2], series.ks_x, [x[2] for x in series.xs]; marker = :circle, label = L"\dot{\theta}", ylabel = L"\dot{\theta}")
    plot!(fig_nom[3], series.ks_u, [u[1] for u in series.us]; marker = :circle, label = L"u", xlabel = "k", ylabel = "u")
    nominal_path = _savefig_pdf(fig_nom, plots_dir, "trajectory_std_nominal_state_control_vs_step.pdf")

    umin, umax = _input_bounds(run_result.problem, cfg)
    uvals = [u[1] for u in series.us]
    nominal_margins = [min(u - umin, umax - u) for u in uvals]
    feedback_margins = Float64[]
    feedback_ks = Int[]
    for step in sort(cert.steps; by = s -> s.k)
        step.status == :ok || continue
        req = _step_field(step, :required_U_radius)
        req === nothing && continue
        uc_margin = step.k <= length(uvals) ? min(uvals[step.k] - umin, umax - uvals[step.k]) : NaN
        push!(feedback_ks, step.k)
        push!(feedback_margins, uc_margin - maximum(Float64.(collect(req))))
    end
    fig_input = plot(;
        xlabel = "k",
        ylabel = "input margin",
        thesis_plot_kwargs(; legend = :topright, size = (700, 420))...,
    )
    plot!(fig_input, series.ks_u, nominal_margins; marker = :circle, label = "nominal input margin")
    plot!(fig_input, feedback_ks, feedback_margins; marker = :diamond, label = "feedback image margin estimate")
    hline!(fig_input, [0.0]; color = :black, ls = :dash, label = "input bound")
    input_margin_path = _savefig_pdf(fig_input, plots_dir, "trajectory_std_input_margin_vs_step.pdf")

    zoom_initial_path = _savefig_pdf(
        _build_zoom_phase_plot(run_result, initial_k; title = "initial ellipsoid window"),
        plots_dir,
        "trajectory_std_zoom_initial.pdf",
    )
    zoom_min_path = _savefig_pdf(
        _build_zoom_phase_plot(run_result, min_k; title = "minimum-volume ellipsoid window"),
        plots_dir,
        "trajectory_std_zoom_min_volume.pdf",
    )
    zoom_terminal_path = _savefig_pdf(
        _build_zoom_phase_plot(run_result, terminal_k; title = "terminal ellipsoid window"),
        plots_dir,
        "trajectory_std_zoom_terminal.pdf",
    )
    selected = unique([initial_k, min_k, clamp(round(Int, 0.33 * terminal_k), initial_k, terminal_k), clamp(round(Int, 0.66 * terminal_k), initial_k, terminal_k)])
    while length(selected) < 4
        push!(selected, terminal_k)
    end
    zoom_selected_path = _savefig_pdf(
        _build_selected_steps_zoom_plot(run_result, selected[1:4]),
        plots_dir,
        "trajectory_std_zoom_selected_steps.pdf",
    )

    println("=== Certification diagnostic summary ===")
    println("smallest volume step k: ", min_k)
    println("smallest volume: ", volumes[min_idx])
    println("minimum X containment margin: ", isempty(xm.vals) ? missing : minimum(xm.vals))
    println("minimum U containment margin: ", isempty(um.vals) ? missing : minimum(um.vals))
    if !isempty(Ls.vals)
        maxL = maximum(maximum(v) for v in Ls.vals)
        println("maximum adaptive remainder L component: ", maxL)
    end
    println("minimum nominal input margin: ", isempty(nominal_margins) ? missing : minimum(nominal_margins))
    println("minimum feedback image input margin estimate: ", isempty(feedback_margins) ? missing : minimum(feedback_margins))

    return (;
        volume_path,
        principal_radii_path,
        remainder_path,
        margins_path,
        box_radii_path,
        nominal_path,
        input_margin_path,
        zoom_initial_path,
        zoom_min_path,
        zoom_terminal_path,
        zoom_selected_path,
        min_volume_step = min_k,
        min_volume = volumes[min_idx],
    )
end

function _stat_summary_or_missing(stat_result, n_stat_samples::Int)
    stat_result === nothing && return (;
        n_samples = n_stat_samples,
        n_reach_target = missing,
        reach_target_rate = missing,
        n_closed_loop_success = missing,
        closed_loop_success_rate = missing,
        n_certified_success = missing,
        certified_success_rate = missing,
        n_violate_safety_or_input = missing,
        violation_rate = missing,
    )
    summary = stat_result.summary
    n_left = hasproperty(summary, :n_left_domain) ? summary.n_left_domain : 0
    n_obs = hasproperty(summary, :n_hit_obstacle) ? summary.n_hit_obstacle : 0
    n_violate = n_left + n_obs
    return (;
        n_samples = summary.n_samples,
        n_reach_target = hasproperty(summary, :n_final_in_target) ? summary.n_final_in_target : missing,
        reach_target_rate = hasproperty(summary, :final_target_rate) ? summary.final_target_rate : missing,
        n_closed_loop_success = hasproperty(summary, :n_closed_loop_success) ?
                                summary.n_closed_loop_success : missing,
        closed_loop_success_rate = hasproperty(summary, :closed_loop_success_rate) ?
                                   summary.closed_loop_success_rate : missing,
        n_certified_success = hasproperty(summary, :n_certified_success) ?
                              summary.n_certified_success : missing,
        certified_success_rate = hasproperty(summary, :certified_success_rate) ?
                                 summary.certified_success_rate : missing,
        n_violate_safety_or_input = n_violate,
        violation_rate = summary.n_samples > 0 ? n_violate / summary.n_samples : missing,
    )
end

function _scaling_string(scaling)
    scaling === nothing && return "nothing"
    return "[" * join(string.(Float64.(scaling)), "; ") * "]"
end

function _display_value(x; digits = 4)
    x === missing && return "missing"
    x === nothing && return "nothing"
    x isa AbstractFloat && return Printf.@sprintf("%.*g", digits, x)
    return string(x)
end

function _csv_value(x)
    text = x === missing || x === nothing ? "" : string(x)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function _rank_key(summary)
    cert_ok = summary.certification_success ? 1 : 0
    formal_ok = summary.formal_box_containment_pass === true ? 1 : 0
    stat_ok = summary.statistical_all_reach_target === true ? 1 : 0
    initial = summary.initial_volume === missing ? -Inf : Float64(summary.initial_volume)
    minvol = summary.minimum_volume === missing ? -Inf : Float64(summary.minimum_volume)
    medvol = summary.median_volume === missing ? -Inf : Float64(summary.median_volume)
    return (cert_ok, formal_ok, stat_ok, initial, minvol, medvol)
end

function _summarize_scaling_run(name, scaling, run_result, stat_result, n_stat_samples)
    cert = run_result.result.certification
    volumes = _certification_volumes(cert)
    box = _adaptive_box_margins(cert)
    min_x_margin = isempty(box.x_margins) ? missing : minimum(box.x_margins)
    min_u_margin = isempty(box.u_margins) ? missing : minimum(box.u_margins)
    formal_pass =
        cert.success &&
        min_x_margin !== missing &&
        min_u_margin !== missing &&
        min_x_margin >= 1.0 - 1.0e-8 &&
        min_u_margin >= 1.0 - 1.0e-8
    stat = _stat_summary_or_missing(stat_result, n_stat_samples)
    statistical_all_reach_target =
        stat.n_reach_target === missing ? missing : stat.n_reach_target == stat.n_samples
    return (;
        candidate_name = name,
        scaling_vector = _scaling_string(scaling),
        certification_success = cert.success,
        certified_transitions = length(cert.steps),
        failed_step = cert.failed_k,
        solve_time = cert.solve_time_sec,
        initial_volume = cert.success && !isempty(volumes) ? last(volumes) : missing,
        minimum_volume = isempty(volumes) ? missing : minimum(volumes),
        median_volume = isempty(volumes) ? missing : Statistics.median(volumes),
        terminal_volume = isempty(volumes) ? missing : first(volumes),
        average_adaptive_iterations = _safe_mean(box.iters),
        maximum_adaptive_iterations = _safe_maximum(box.iters),
        average_candidate_boxes = _safe_mean(box.candidates),
        maximum_candidate_boxes = _safe_maximum(box.candidates),
        min_X_containment_margin = min_x_margin,
        min_U_containment_margin = min_u_margin,
        formal_box_containment_pass = formal_pass,
        n_statistical_samples = stat.n_samples,
        n_reach_target = stat.n_reach_target,
        reach_target_rate = stat.reach_target_rate,
        n_closed_loop_success = stat.n_closed_loop_success,
        closed_loop_success_rate = stat.closed_loop_success_rate,
        n_certified_success = stat.n_certified_success,
        certified_success_rate = stat.certified_success_rate,
        n_violate_safety_or_input = stat.n_violate_safety_or_input,
        violation_rate = stat.violation_rate,
        statistical_all_reach_target = statistical_all_reach_target,
    )
end

function _print_scaling_sweep_table(rows, best)
    println("\n=== Simple pendulum scaling sweep: adaptive boxes, max-volume ===")
    println("Scaling alone is not accepted unless post-hoc containment checks pass.")
    header = Printf.@sprintf(
        "%-20s %-6s %-6s %-6s %12s %12s %12s %10s %10s %8s",
        "candidate",
        "cert",
        "boxes",
        "stats",
        "V0",
        "Vmin",
        "Vmed",
        "Xmargin",
        "Umargin",
        "time",
    )
    println(header)
    println("-" ^ length(header))
    for row in rows
        stats_ok = row.statistical_all_reach_target === true ? "yes" :
                   row.statistical_all_reach_target === false ? "no" : "miss"
        println(
            Printf.@sprintf(
                "%-20s %-6s %-6s %-6s %12s %12s %12s %10s %10s %8s",
                row.candidate_name,
                row.certification_success ? "yes" : "no",
                row.formal_box_containment_pass ? "yes" : "no",
                stats_ok,
                _display_value(row.initial_volume),
                _display_value(row.minimum_volume),
                _display_value(row.median_volume),
                _display_value(row.min_X_containment_margin),
                _display_value(row.min_U_containment_margin),
                _display_value(row.solve_time),
            ),
        )
    end
    println("\nBest valid configuration: ", best === nothing ? "none" : best.candidate_name)
    return nothing
end

function _write_scaling_sweep_csv(path, rows)
    mkpath(dirname(path))
    columns = [
        :candidate_name,
        :scaling_vector,
        :certification_success,
        :certified_transitions,
        :failed_step,
        :solve_time,
        :initial_volume,
        :minimum_volume,
        :median_volume,
        :terminal_volume,
        :average_adaptive_iterations,
        :maximum_adaptive_iterations,
        :average_candidate_boxes,
        :maximum_candidate_boxes,
        :min_X_containment_margin,
        :min_U_containment_margin,
        :formal_box_containment_pass,
        :n_statistical_samples,
        :n_reach_target,
        :reach_target_rate,
        :n_closed_loop_success,
        :closed_loop_success_rate,
        :n_certified_success,
        :certified_success_rate,
        :n_violate_safety_or_input,
        :violation_rate,
        :statistical_all_reach_target,
    ]
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            println(io, join((_csv_value(getproperty(row, col)) for col in columns), ","))
        end
    end
    return path
end

function _write_scaling_sweep_markdown(path, rows, best)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Simple pendulum scaling sweep")
        println(io)
        println(io, "Adaptive linearization boxes were enabled with objective `:max_volume`.")
        println(io, "Scaling alone is not accepted unless post-hoc containment checks pass.")
        println(io)
        println(io, "Best valid configuration: ", best === nothing ? "none" : best.candidate_name)
        println(io)
        println(
            io,
            "| candidate | cert | boxes | all target | V0 | Vmin | Vmed | X margin | U margin | samples |",
        )
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in rows
            println(
                io,
                "| ",
                row.candidate_name,
                " | ",
                row.certification_success,
                " | ",
                row.formal_box_containment_pass,
                " | ",
                row.statistical_all_reach_target,
                " | ",
                _display_value(row.initial_volume),
                " | ",
                _display_value(row.minimum_volume),
                " | ",
                _display_value(row.median_volume),
                " | ",
                _display_value(row.min_X_containment_margin),
                " | ",
                _display_value(row.min_U_containment_margin),
                " | ",
                row.n_statistical_samples,
                " |",
            )
        end
    end
    return path
end

function print_adaptive_box_step_diagnostics(cert)
    println("=== Linearization box step diagnostics ===")
    for step in cert.steps
        step.status == :ok || continue
        println(
            "k=",
            step.k,
            " status=",
            _step_field(step, :adaptive_box_status),
            " selected_scale=",
            _step_field(step, :selected_scale),
            " selected_logvolume=",
            _step_field(step, :selected_logvolume),
            " δx=",
            _step_field(step, :Xbar_radius),
            " δu=",
            _step_field(step, :Ubar_radius),
            " required_δx=",
            _step_field(step, :required_X_radius),
            " required_δu=",
            _step_field(step, :required_U_radius),
            " L=",
            _step_field(step, :L),
            " adaptive_iters=",
            _step_field(step, :adaptive_box_iters),
            " candidate_boxes=",
            _step_field(step, :number_of_candidate_boxes),
        )
    end
    return nothing
end

function print_mppi_terminal_diagnostics(run_result, cfg, terminal_data)
    println("=== MPPI terminal diagnostics ===")
    xs = collect(ST.enum_elems(run_result.nominal_candidate.x_traj))
    xN = Float64.(collect(xs[end]))
    println("nominal endpoint x_N: ", xN)
    if terminal_data === nothing
        println("terminal ellipsoidal distance dT2(x_N): NaN")
        println("nominal endpoint inside terminal ellipsoid: false")
        return nothing
    end
    dT2 = terminal_ellipsoidal_distance2(
        xN,
        terminal_data.terminal_center,
        terminal_data.terminal_shape,
    )
    println("terminal ellipsoidal distance dT2(x_N): ", dT2)
    println("nominal endpoint inside terminal ellipsoid: ", dT2 <= 1.0 + 1.0e-8)
    return nothing
end

function print_certification_diagnostics(run_result, cfg, terminal_data)
    cert = run_result.result.certification
    println("=== Certification diagnostics ===")
    println("certification success: ", cert.success)
    println("number of certified transitions: ", length(cert.steps))
    println("failed step: ", cert.failed_k)
    println("runtime: ", cert.solve_time_sec)
    println(
        "terminal ellipsoid used by certifier: ",
        terminal_data === nothing ? "radius fallback" : "fixed inner ellipsoid",
    )
    println("state scaling used by certifier: ", certifier_state_scaling(cfg))
    println("adaptive linearization boxes: ", cfg.adaptive_linearization_boxes)
    println("adaptive box objective: ", cfg.adaptive_box_objective)
    println("adaptive box keep first consistent: ", cfg.adaptive_box_keep_first_consistent)
    println("adaptive box search scales: ", cfg.adaptive_box_search_scales)
    print_adaptive_box_step_diagnostics(cert)
    return nothing
end

function main(
    cfg::SimplePendulumMPPIConfig = SimplePendulumMPPIConfig();
    print_diagnostics::Bool = true,
    save_outputs::Bool = true,
    artifact_prefix::AbstractString = "",
)
    Random.seed!(cfg.rng_seed)

    certification_basename =
        isempty(artifact_prefix) ? "ellipsoids" : "$(artifact_prefix)_certification"
    warmup_basename =
        isempty(artifact_prefix) ? "abstract_traj_used_as_warmup" : "$(artifact_prefix)_warmup"
    nominal_basename =
        isempty(artifact_prefix) ? "mppi_candidate_traj" : "$(artifact_prefix)_nominal"
    statistical_basename =
        isempty(artifact_prefix) ? "kappa_statistical_rollouts" : "$(artifact_prefix)_statistical"

    input_mapping = build_pendulum_input_mapping(cfg)

    generator_builder = (problem, system_cfg, control_cfg, cfg_) ->
        build_simple_pendulum_mppi_generator(
            problem,
            system_cfg,
            control_cfg,
            cfg_;
            input_mapping = input_mapping,
        )

    certifier_builder = (problem, system_cfg, control_cfg, cfg_) ->
        build_simple_pendulum_certifier(problem, system_cfg, control_cfg, cfg_)

    save_artifacts! = function (run_result)
        save_outputs || return (;)
        save_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = certification_basename,
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            phase_title = "simple_pendulum_mppi - certification phase space",
            state_title = "simple_pendulum_mppi - certification states",
            control_title = "simple_pendulum_mppi - certification control",
        )
        animation_paths = save_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = certification_basename,
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "simple_pendulum_mppi - certification",
        )
        return (; gif = animation_paths.gif_path, mp4 = animation_paths.mp4_path)
    end

    run_result = run_benchmark(
        cfg;
        scenario_name = "simple_pendulum_mppi",
        build_concrete_system = () -> build_pendulum_system_cfg(cfg; pendulum_module = SP),
        build_control_problem = () -> build_pendulum_control_cfg(cfg; pendulum_module = SP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )

    terminal_data = terminal_inner_ellipsoid_data(run_result.problem, cfg)
    if print_diagnostics
        print_terminal_inner_ellipsoid_diagnostics(run_result.problem, cfg, terminal_data)
        print_mppi_terminal_diagnostics(run_result, cfg, terminal_data)
        print_certification_diagnostics(run_result, cfg, terminal_data)
    end

    warmup_candidate = OP.get_seed(run_result.solver.generator)
    warmup_candidate === nothing && error("MPPI warmup trajectory is missing.")

    if save_outputs
        save_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            warmup_candidate;
            basename = warmup_basename,
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            phase_title = "simple_pendulum_mppi - warmup phase space",
            state_title = "simple_pendulum_mppi - warmup states",
            control_title = "simple_pendulum_mppi - warmup control",
        )
        save_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            warmup_candidate;
            basename = warmup_basename,
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "simple_pendulum_mppi - warmup",
        )

        save_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = nominal_basename,
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            phase_title = "simple_pendulum_mppi - candidate phase space",
            state_title = "simple_pendulum_mppi - candidate states",
            control_title = "simple_pendulum_mppi - candidate control",
        )
        save_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            run_result.nominal_candidate;
            basename = nominal_basename,
            cert_result = nothing,
            show_ellipsoids = false,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "simple_pendulum_mppi - candidate",
        )
    end

    stat_result = nothing
    if run_result.result.certification.success && cfg.kappa_statistical_samples > 0
        stat_result = run_kappa_statistical_check(
            run_result;
            n_samples = cfg.kappa_statistical_samples,
        )
        if save_outputs
            save_kappa_statistical_plots!(
                stat_result;
                basename = statistical_basename,
                wrap_angles = true,
                show_ellipsoids = false,
                axis_labels_12 = (L"\theta\,[\mathrm{rad}]", L"\dot{\theta}\,[\mathrm{rad/s}]"),
                rollout_alpha = 0.10,
                rollout_lw = 0.35,
                rollout_ms = 0.7,
                initial_sample_alpha = 0.22,
                initial_sample_ms = 0.8,
            )
        end
    elseif run_result.result.certification.success
        println("Skipping kappa statistical check because kappa_statistical_samples <= 0.")
    else
        println("Skipping kappa statistical check because certification failed.")
    end

    return merge(run_result, (; statistical_result = stat_result))
end

function run_simple_pendulum_certification_mode(; adaptive_linearization_boxes::Bool)
    cfg = SimplePendulumMPPIConfig(;
        adaptive_linearization_boxes = adaptive_linearization_boxes,
        plot_gif = false,
        plot_mp4 = false,
        kappa_statistical_samples = 0,
    )
    run_result = main(cfg)
    println(
        "mode adaptive_linearization_boxes=",
        adaptive_linearization_boxes,
        " success=",
        run_result.result.certification.success,
    )
    return run_result
end

function compare_fixed_and_adaptive_simple_pendulum()
    fixed = run_simple_pendulum_certification_mode(; adaptive_linearization_boxes = false)
    adaptive = run_simple_pendulum_certification_mode(; adaptive_linearization_boxes = true)
    return (; fixed, adaptive)
end

function run_simple_pendulum_trajectory_std(;
    base_cfg::SimplePendulumMPPIConfig = SimplePendulumMPPIConfig(),
    kappa_statistical_samples::Int = max(base_cfg.kappa_statistical_samples, 12000),
)
    println("=== Simple pendulum trajectory_std final run ===")
    println("First pass: generate nominal MPPI trajectory and compute trajectory_std scaling.")

    first_pass_cfg = _with_config(
        base_cfg;
        use_state_scaling = false,
        adaptive_linearization_boxes = true,
        adaptive_box_objective = :max_volume,
        adaptive_box_keep_first_consistent = false,
        adaptive_box_verbose = false,
        plot_gif = false,
        plot_mp4 = false,
        kappa_statistical_samples = 0,
    )
    first_pass = main(first_pass_cfg; print_diagnostics = false, save_outputs = false)

    scaling = compute_trajectory_std_scaling_from_candidate(
        first_pass.certification_candidate,
        first_pass_cfg,
    )
    println("computed trajectory_std state scaling: ", scaling)
    println("Second pass: final certification with trajectory_std scaling and adaptive boxes.")

    final_cfg = _with_config(
        base_cfg;
        use_state_scaling = true,
        state_scaling = Float64.(scaling),
        adaptive_linearization_boxes = true,
        adaptive_box_objective = :max_volume,
        adaptive_box_keep_first_consistent = false,
        adaptive_box_verbose = false,
        plot_gif = true,
        plot_mp4 = true,
        kappa_statistical_samples = kappa_statistical_samples,
    )
    final_run = main(
        final_cfg;
        print_diagnostics = true,
        save_outputs = true,
        artifact_prefix = "trajectory_std",
    )

    cert = final_run.result.certification
    volumes = _certification_volume_summary(cert)
    n_ellipsoids = length(extract_ellipsoids(cert; max_keep = nothing))
    diagnostics = save_simple_pendulum_certification_diagnostics!(
        final_run,
        final_cfg;
        basename = "trajectory_std_diagnostics",
    )
    println("=== Trajectory_std final summary ===")
    println("computed trajectory_std scaling: ", scaling)
    println("certification success: ", cert.success)
    println("number of certified transitions: ", length(cert.steps))
    println("certified ellipsoids plotted: ", n_ellipsoids)
    println("initial ellipsoid volume: ", volumes.initial)
    println("minimum ellipsoid volume: ", volumes.minimum)
    println("median ellipsoid volume: ", volumes.median)
    println("plots directory: ", final_run.outputs.plots_dir)
    println("animations directory: ", final_run.outputs.animations_dir)
    println("All certification plots use every ellipsoid extracted from the certification chain.")

    return (;
        scaling,
        first_pass,
        final_run,
        volume_summary = volumes,
        n_ellipsoids,
        diagnostics,
    )
end

function _run_trajectory_std_certification_for_cfg(
    base_cfg::SimplePendulumMPPIConfig;
    n_stat_samples::Int = 0,
    save_outputs::Bool = false,
    print_diagnostics::Bool = false,
    artifact_prefix::AbstractString = "",
)
    first_cfg = _with_config(
        base_cfg;
        use_state_scaling = false,
        adaptive_linearization_boxes = true,
        adaptive_box_objective = :max_volume,
        adaptive_box_keep_first_consistent = false,
        adaptive_box_verbose = false,
        plot_gif = false,
        plot_mp4 = false,
        kappa_statistical_samples = 0,
    )
    first = main(first_cfg; print_diagnostics = false, save_outputs = false)
    scaling = compute_trajectory_std_scaling_from_candidate(first.certification_candidate, first_cfg)
    final_cfg = _with_config(
        base_cfg;
        use_state_scaling = true,
        state_scaling = Float64.(scaling),
        adaptive_linearization_boxes = true,
        adaptive_box_objective = :max_volume,
        adaptive_box_keep_first_consistent = false,
        adaptive_box_verbose = false,
        plot_gif = save_outputs && base_cfg.plot_gif,
        plot_mp4 = save_outputs && base_cfg.plot_mp4,
        kappa_statistical_samples = n_stat_samples,
    )
    final = main(
        final_cfg;
        print_diagnostics = print_diagnostics,
        save_outputs = save_outputs,
        artifact_prefix = artifact_prefix,
    )
    return (; first, final, final_cfg, scaling)
end

function _trajectory_std_summary_row(
    name::AbstractString,
    scaling,
    run_result,
    stat_result;
    extra = (;),
    n_stat_samples::Int = 0,
)
    cert = run_result.result.certification
    volumes = _certification_volume_summary(cert)
    box = _adaptive_box_margins(cert)
    min_x_margin = isempty(box.x_margins) ? missing : minimum(box.x_margins)
    min_u_margin = isempty(box.u_margins) ? missing : minimum(box.u_margins)
    formal_pass =
        cert.success &&
        min_x_margin !== missing &&
        min_u_margin !== missing &&
        min_x_margin >= 1.0 - 1.0e-8 &&
        min_u_margin >= 1.0 - 1.0e-8
    stat = _stat_summary_or_missing(stat_result, n_stat_samples)
    statistical_all_reach_target =
        stat.n_reach_target === missing ? missing : stat.n_reach_target == stat.n_samples
    row = (;
        candidate_name = name,
        scaling_vector = _scaling_string(scaling),
        certification_success = cert.success,
        certified_transitions = length(cert.steps),
        failed_step = cert.failed_k,
        solve_time = cert.solve_time_sec,
        initial_volume = volumes.initial,
        minimum_volume = volumes.minimum,
        median_volume = volumes.median,
        terminal_volume = volumes.terminal,
        min_X_containment_margin = min_x_margin,
        min_U_containment_margin = min_u_margin,
        formal_box_containment_pass = formal_pass,
        n_statistical_samples = stat.n_samples,
        n_reach_target = stat.n_reach_target,
        reach_target_rate = stat.reach_target_rate,
        n_closed_loop_success = stat.n_closed_loop_success,
        closed_loop_success_rate = stat.closed_loop_success_rate,
        n_violate_safety_or_input = stat.n_violate_safety_or_input,
        violation_rate = stat.violation_rate,
        statistical_all_reach_target = statistical_all_reach_target,
    )
    return merge(extra, row)
end

function _generic_rank_key(row)
    cert_ok = row.certification_success ? 1 : 0
    formal_ok = row.formal_box_containment_pass === true ? 1 : 0
    stat_ok = row.statistical_all_reach_target === true ? 1 : 0
    minvol = row.minimum_volume === missing ? -Inf : Float64(row.minimum_volume)
    initial = row.initial_volume === missing ? -Inf : Float64(row.initial_volume)
    medvol = row.median_volume === missing ? -Inf : Float64(row.median_volume)
    return (cert_ok, formal_ok, stat_ok, minvol, initial, medvol)
end

function _write_rows_csv(path, rows)
    mkpath(dirname(path))
    columns = Symbol[]
    for row in rows
        for name in propertynames(row)
            name in columns || push!(columns, name)
        end
    end
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            println(io, join((_csv_value(hasproperty(row, col) ? getproperty(row, col) : missing) for col in columns), ","))
        end
    end
    return path
end

function _write_rows_report(path, title, rows, best)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        println(io, "Best valid configuration: ", best === nothing ? "none" : best.candidate_name)
        println(io)
        println(io, "| candidate | cert | boxes | all target | Vmin | V0 | Vmed | X margin | U margin |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in rows
            println(
                io,
                "| ", row.candidate_name,
                " | ", row.certification_success,
                " | ", row.formal_box_containment_pass,
                " | ", row.statistical_all_reach_target,
                " | ", _display_value(row.minimum_volume),
                " | ", _display_value(row.initial_volume),
                " | ", _display_value(row.median_volume),
                " | ", _display_value(row.min_X_containment_margin),
                " | ", _display_value(row.min_U_containment_margin),
                " |",
            )
        end
    end
    return path
end

function _best_valid_row(rows)
    ranked = sort(rows; by = _generic_rank_key, rev = true)
    for row in ranked
        if row.certification_success &&
           row.formal_box_containment_pass &&
           row.statistical_all_reach_target === true
            return ranked, row
        end
    end
    return ranked, nothing
end

function sweep_simple_pendulum_terminal_shrink_trajectory_std(;
    base_cfg::SimplePendulumMPPIConfig = SimplePendulumMPPIConfig(),
    terminal_shrinks = [0.85, 0.90, 0.95, 1.00],
    n_stat_samples::Int = 300,
)
    rows = NamedTuple[]
    for shrink in terminal_shrinks
        println("\n--- terminal_shrink trajectory_std sweep: ", shrink, " ---")
        cfg = _with_config(
            base_cfg;
            terminal_shrink = Float64(shrink),
            plot_gif = false,
            plot_mp4 = false,
        )
        result = _run_trajectory_std_certification_for_cfg(
            cfg;
            n_stat_samples = n_stat_samples,
            save_outputs = false,
            print_diagnostics = false,
        )
        row = _trajectory_std_summary_row(
            "terminal_shrink=$(shrink)",
            result.scaling,
            result.final,
            result.final.statistical_result;
            extra = (; terminal_shrink = Float64(shrink)),
            n_stat_samples = n_stat_samples,
        )
        push!(rows, row)
    end

    ranked, best = _best_valid_row(rows)
    csv_path = joinpath(base_cfg.output_root, "simple_pendulum_terminal_shrink_sweep_results.csv")
    md_path = joinpath(base_cfg.output_root, "simple_pendulum_terminal_shrink_sweep_report.md")
    _write_rows_csv(csv_path, ranked)
    _write_rows_report(md_path, "Simple pendulum terminal_shrink trajectory_std sweep", ranked, best)
    println("terminal_shrink sweep CSV report written to: ", csv_path)
    println("terminal_shrink sweep Markdown report written to: ", md_path)
    return (; rows = ranked, best, csv_path, markdown_path = md_path)
end

function sweep_simple_pendulum_funnel_friendly_mppi(;
    base_cfg::SimplePendulumMPPIConfig = SimplePendulumMPPIConfig(),
    nsteps = [55, 70, 90],
    velocity_running_weights = [0.15, 0.3, 0.6],
    input_margin_weights = [0.0, 1.0, 5.0],
    terminal_center_weights = [1.0e5, 3.0e5, 1.0e6],
    n_stat_samples::Int = 300,
)
    rows = NamedTuple[]
    for nstep in nsteps,
        velocity_weight in velocity_running_weights,
        input_weight in input_margin_weights,
        terminal_weight in terminal_center_weights

        name = "n=$(nstep)_v=$(velocity_weight)_um=$(input_weight)_tc=$(terminal_weight)"
        println("\n--- funnel-friendly MPPI sweep: ", name, " ---")
        cfg = _with_config(
            base_cfg;
            nstep = Int(nstep),
            velocity_running_weight = Float64(velocity_weight),
            input_margin_weight = Float64(input_weight),
            terminal_center_weight = Float64(terminal_weight),
            plot_gif = false,
            plot_mp4 = false,
        )
        result = _run_trajectory_std_certification_for_cfg(
            cfg;
            n_stat_samples = n_stat_samples,
            save_outputs = false,
            print_diagnostics = false,
        )
        row = _trajectory_std_summary_row(
            name,
            result.scaling,
            result.final,
            result.final.statistical_result;
            extra = (;
                nstep = Int(nstep),
                velocity_running_weight = Float64(velocity_weight),
                input_margin_weight = Float64(input_weight),
                terminal_center_weight = Float64(terminal_weight),
            ),
            n_stat_samples = n_stat_samples,
        )
        push!(rows, row)
    end

    ranked, best = _best_valid_row(rows)
    csv_path = joinpath(base_cfg.output_root, "simple_pendulum_funnel_friendly_mppi_sweep_results.csv")
    md_path = joinpath(base_cfg.output_root, "simple_pendulum_funnel_friendly_mppi_sweep_report.md")
    _write_rows_csv(csv_path, ranked)
    _write_rows_report(md_path, "Simple pendulum funnel-friendly MPPI sweep", ranked, best)
    println("funnel-friendly MPPI sweep CSV report written to: ", csv_path)
    println("funnel-friendly MPPI sweep Markdown report written to: ", md_path)
    return (; rows = ranked, best, csv_path, markdown_path = md_path)
end

function sweep_simple_pendulum_scaling_adaptive_boxes(;
    n_stat_samples::Int = 1000,
    base_cfg::SimplePendulumMPPIConfig = SimplePendulumMPPIConfig(),
)
    n_stat_samples >= 0 || throw(ArgumentError("n_stat_samples must be nonnegative"))

    sweep_base_cfg = _with_config(
        base_cfg;
        adaptive_linearization_boxes = true,
        adaptive_box_objective = :max_volume,
        adaptive_box_keep_first_consistent = false,
        adaptive_box_verbose = false,
        plot_gif = false,
        plot_mp4 = false,
        kappa_statistical_samples = 0,
    )

    println("=== Simple pendulum scaling sweep setup ===")
    println("adaptive_linearization_boxes = true")
    println("adaptive_box_objective = :max_volume")
    println("adaptive_box_keep_first_consistent = false")
    println("Scaling alone is not accepted unless post-hoc containment checks pass.")

    baseline_cfg = _with_config(
        sweep_base_cfg;
        use_state_scaling = false,
        state_scaling = Float64.(sweep_base_cfg.state_scaling),
    )
    println("\n--- scaling candidate: no_scaling ---")
    baseline_run = main(baseline_cfg; print_diagnostics = false)

    candidates = scaling_candidates_from_nominal(baseline_run, sweep_base_cfg)
    rows = NamedTuple[]
    run_results = Dict{String, Any}()
    stat_results = Dict{String, Any}()

    for candidate in candidates
        cfg = _with_config(
            sweep_base_cfg;
            use_state_scaling = candidate.use_state_scaling,
            state_scaling = candidate.scaling === nothing ?
                            Float64.(sweep_base_cfg.state_scaling) :
                            Float64.(candidate.scaling),
            rng_seed = sweep_base_cfg.rng_seed,
            kappa_statistical_samples = 0,
        )

        run_result = if candidate.name == "no_scaling"
            baseline_run
        else
            println("\n--- scaling candidate: ", candidate.name, " ---")
            println("state_scaling = ", candidate.scaling)
            main(cfg; print_diagnostics = false)
        end
        run_results[candidate.name] = run_result

        stat_result = nothing
        if run_result.result.certification.success && n_stat_samples > 0
            stat_result = run_kappa_statistical_check(
                run_result;
                n_samples = n_stat_samples,
                rng = Random.MersenneTwister(sweep_base_cfg.rng_seed + 100_000),
                verbose = true,
            )
        elseif !run_result.result.certification.success
            println("Skipping statistical check for ", candidate.name, " because certification failed.")
        else
            println("Skipping statistical check for ", candidate.name, " because n_stat_samples == 0.")
        end
        stat_results[candidate.name] = stat_result

        row = _summarize_scaling_run(
            candidate.name,
            candidate.scaling,
            run_result,
            stat_result,
            n_stat_samples,
        )
        push!(rows, row)
    end

    ranked = sort(rows; by = _rank_key, rev = true)
    best = nothing
    for row in ranked
        if row.certification_success &&
           row.formal_box_containment_pass &&
           row.statistical_all_reach_target === true
            best = row
            break
        end
    end

    _print_scaling_sweep_table(ranked, best)

    output_dir = joinpath(sweep_base_cfg.output_root)
    csv_path = joinpath(output_dir, "simple_pendulum_scaling_sweep_results.csv")
    md_path = joinpath(output_dir, "simple_pendulum_scaling_sweep_report.md")
    _write_scaling_sweep_csv(csv_path, ranked)
    _write_scaling_sweep_markdown(md_path, ranked, best)
    println("Scaling sweep CSV report written to: ", csv_path)
    println("Scaling sweep Markdown report written to: ", md_path)

    return (;
        rows = ranked,
        best = best,
        csv_path = csv_path,
        markdown_path = md_path,
        run_results = run_results,
        stat_results = stat_results,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_simple_pendulum_trajectory_std()
end
