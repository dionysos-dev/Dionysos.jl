include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
using LaTeXStrings
import LinearAlgebra as LA
import MathematicalSystems as MS
import Random

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "pendulum",
        "simple_pendulum.jl",
    ),
)
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
    nstep::Int = 75
    input_values::Tuple{Vararg{Float64}} = Tuple(-2.3:0.25:2.3)

    terminal_radius::Float64 = 0.25
    use_terminal_inner_ellipsoid::Bool = true
    terminal_shrink::Float64 = 0.85
    terminal_outside_weight::Float64 = 1.0e6
    terminal_center_weight::Float64 = 1.0e5
    terminal_success_distance2::Float64 = 0.1
    λ::Float64 = 0.005
    maxδx::Float64 = 80.0
    maxδu::Float64 = 40.0
    use_state_scaling::Bool = true
    state_scaling::Vector{Float64} = [0.85 * 15.0 * pi / 180.0, 0.25]
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-0.02, 0.02), IA.interval(-0.02, 0.02))
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.02, 0.02), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)
    adaptive_linearization_boxes::Bool = true
    ΔX_initial::Vector{Float64} = [0.02, 0.02]
    ΔX_min::Vector{Float64} = [1.0e-8, 1.0e-8]
    ΔX_max::Vector{Float64} = [0.30, 2.00]
    ΔU_initial::Vector{Float64} = [0.02]
    ΔU_min::Vector{Float64} = [1.0e-3]
    ΔU_max::Vector{Float64} = [1.00]
    adaptive_box_growth::Float64 = 1.5
    adaptive_box_safety::Float64 = 1.01
    adaptive_box_max_iters::Int = 16
    adaptive_box_atol::Float64 = 1.0e-8
    adaptive_box_verbose::Bool = true
    adaptive_box_search_scales::Vector{Float64} =
        [0.5, 0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 2.0]
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

    planning_input_scale::Float64 = 1.0
    mppi_nsamples::Int = 7800
    mppi_niter::Int = 20
    mppi_λ::Float64 = 1.75
    mppi_noise_u::Float64 = 1.8
    rng_seed::Int = 1
    kappa_statistical_samples::Int = 300
end

function build_inner_terminal_ellipsoid(lower, upper; shrink = 0.95)
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
        throw(ArgumentError("terminal_shrink must be finite and in (0, 1]"))

    radii = 0.5 .* (ub .- lb)
    all(radii .> 0.0) || throw(ArgumentError("target radii must be positive"))

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
    return all(center .+ radii .<= ub .+ atol) && all(center .- radii .>= lb .- atol)
end

terminal_ellipsoidal_distance2(x, terminal_center, terminal_shape) = begin
    e = Float64.(collect(x)) .- terminal_center
    return Float64(e' * terminal_shape * e)
end

function input_bounds(problem, cfg::SimplePendulumMPPIConfig)
    if hasproperty(problem.system, :U) && problem.system.U isa UT.HyperRectangle
        return Float64(problem.system.U.lb[1]), Float64(problem.system.U.ub[1])
    end
    vals = Float64.(collect(cfg.input_values))
    return minimum(vals), maximum(vals)
end

function planning_input_bounds(problem, cfg::SimplePendulumMPPIConfig)
    0.0 < cfg.planning_input_scale <= 1.0 ||
        throw(ArgumentError("planning_input_scale must be in (0, 1]"))
    umin, umax = input_bounds(problem, cfg)
    center = 0.5 * (umin + umax)
    half_width = 0.5 * (umax - umin)
    return center - cfg.planning_input_scale * half_width,
    center + cfg.planning_input_scale * half_width
end

function planning_input_domain(problem, cfg::SimplePendulumMPPIConfig)
    umin_plan, umax_plan = planning_input_bounds(problem, cfg)
    return UT.HyperRectangle([umin_plan], [umax_plan])
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
    terminal_center, terminal_shape =
        build_inner_terminal_ellipsoid(problem.target_set; shrink = cfg.terminal_shrink)
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
    # MPPI plans with U_plan, a restricted subset of the physical U_cert.
    # The certifier below still receives problem.system.U, i.e. the full U_cert.
    plan_U = planning_input_domain(problem, cfg)
    project_input = u -> project_pendulum_input_to_domain(u, plan_U)

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
            J += 0.15 * e_goal[2]^2

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
            hit_idx = findfirst(
                d -> d <= cfg.terminal_success_distance2 + 1.0e-8,
                terminal_distances,
            )
            hit_index = hit_idx === nothing ? argmin(terminal_distances) : hit_idx
            hit_target = hit_idx !== nothing
        end

        xT = wrap_state(xs[min(hit_index, length(xs))])
        eT = pendulum_state_error(xT, target_center, cfg)
        J += 80.0 * eT[1]^2
        J += 8.0 * eT[2]^2

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
        for k in 1:last_u_index
            uk = collect(Float64, us[k])
            J += 0.02 * uk[1]^2

            if k >= 2
                duk = uk[1] - us[k - 1][1]
                J += 0.04 * duk^2
            end
        end

        return J
    end

    success_fun = if terminal_data === nothing
        (prob, cand) ->
            any(x -> (wrap_state(x) ∈ prob.target_set), ST.enum_elems(cand.x_traj))
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
        (_prob, _mppi_cfg, cand) -> truncate_candidate_at_first_terminal_ellipsoid_hit(
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
        terminal_center = terminal_data === nothing ? nothing :
                          terminal_data.terminal_center,
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
    save_outputs::Bool = true,
    run_statistical::Bool = true,
    artifact_prefix::AbstractString = "",
)
    Random.seed!(cfg.rng_seed)

    certification_basename =
        isempty(artifact_prefix) ? "ellipsoids" : "$(artifact_prefix)_certification"
    warmup_basename =
        isempty(artifact_prefix) ? "abstract_traj_used_as_warmup" :
        "$(artifact_prefix)_warmup"
    nominal_basename =
        isempty(artifact_prefix) ? "mppi_candidate_traj" : "$(artifact_prefix)_nominal"

    input_mapping = build_pendulum_input_mapping(cfg)

    generator_builder =
        (problem, system_cfg, control_cfg, cfg_) -> build_simple_pendulum_mppi_generator(
            problem,
            system_cfg,
            control_cfg,
            cfg_;
            input_mapping = input_mapping,
        )

    certifier_builder =
        (problem, system_cfg, control_cfg, cfg_) ->
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
        build_concrete_system = () ->
            build_pendulum_system_cfg(cfg; pendulum_module = SP),
        build_control_problem = () ->
            build_pendulum_control_cfg(cfg; pendulum_module = SP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )

    terminal_data = terminal_inner_ellipsoid_data(run_result.problem, cfg)
    print_terminal_inner_ellipsoid_diagnostics(run_result.problem, cfg, terminal_data)
    print_mppi_terminal_diagnostics(run_result, cfg, terminal_data)
    print_certification_diagnostics(run_result, cfg, terminal_data)

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
    if run_statistical &&
       run_result.result.certification.success &&
       cfg.kappa_statistical_samples > 0
        stat_result = run_kappa_statistical_check(
            run_result;
            n_samples = cfg.kappa_statistical_samples,
        )
        if save_outputs
            save_kappa_statistical_plots!(
                stat_result;
                output_dir = run_result.outputs.plots_dir,
                basename = isempty(artifact_prefix) ? "kappa_statistical_rollouts" :
                           "$(artifact_prefix)_statistical",
                wrap_angles = true,
                axis_labels_12 = (
                    L"\theta\,[\mathrm{rad}]",
                    L"\dot{\theta}\,[\mathrm{rad/s}]",
                ),
            )
        end
    elseif run_statistical && run_result.result.certification.success
        println("Skipping kappa statistical check because kappa_statistical_samples <= 0.")
    elseif run_statistical
        println("Skipping kappa statistical check because certification failed.")
    end

    return merge(run_result, (; statistical = stat_result))
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

if abspath(PROGRAM_FILE) == @__FILE__
    compare_fixed_and_adaptive_simple_pendulum()
end
