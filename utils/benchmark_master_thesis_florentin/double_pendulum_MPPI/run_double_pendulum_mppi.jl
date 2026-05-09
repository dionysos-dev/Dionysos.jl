include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import MathematicalSystems as MS
import Random

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pendulum", "double_pendulum.jl"))
const DP = DoublePendulum

Base.@kwdef struct DoublePendulumMPPIConfig
    l1::Float64 = 1.0
    l2::Float64 = 1.0
    m1::Float64 = 1.0
    m2::Float64 = 1.0
    g::Float64 = 9.81
    objective::String = "benchmark_up_convex"

    Δt::Float64 = 0.11
    hx::SVector{4, Float64} = SVector(5 * (pi / 180), 5 * (pi / 180), 0.25, 0.25)
    periodic_dims::SVector{2, Int} = SVector(1, 2)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 41
    input_values::Tuple{Vararg{Float64}} = vals_tuple = Tuple(-3.5:0.25:3.5)

    terminal_radius::Float64 = 0.3
    λ::Float64 = 0.0000005
    maxδx::Float64 = 120.0
    maxδu::Float64 = 60.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
        IA.interval(-0.01, 0.01),
    )
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.01, 0.01), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    plot_mp4::Bool = true
    verbose::Bool = false

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 1

    mppi_nsamples::Int = 1800
    mppi_niter::Int = 20
    mppi_λ::Float64 = 1.75
    mppi_noise_u::Float64 = 0.8
    rng_seed::Int = 3
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

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)
    noise_sampler = (rng, u, k) -> [cfg.mppi_noise_u * Random.randn(rng)]
    project_input = u -> project_input_to_domain(u, problem.system.U)

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

            if xw ∈ prob.target_set
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

        xT = wrap_state(xs[min(hit_index, length(xs))])
        eT = periodic_state_error(xT, target_center, cfg)
        J += 45.0 * eT[1]^2
        J += 90.0 * eT[2]^2
        J += 4.0 * eT[3]^2
        J += 4.0 * eT[4]^2

        if !hit_target
            J += 250.0
        end

        last_u_index = min(length(us), max(hit_index - 1, 0))
        for k in 1:last_u_index
            uk = collect(Float64, us[k])
            J += 0.015 * uk[1]^2

            if k >= 2
                duk = uk[1] - us[k - 1][1]
                J += 0.05 * duk^2
            end
        end

        return J
    end

    success_fun =
        (prob, cand) -> any(x -> (wrap_state(x) ∈ prob.target_set), ST.enum_elems(cand.x_traj))

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
        (prob, x) -> wrap_state(x),
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function build_double_pendulum_certifier(
    problem,
    _system_cfg,
    _control_cfg,
    cfg::DoublePendulumMPPIConfig,
)
    symbolic_builder = build_double_pendulum_symbolic_builder(cfg; pendulum_module = DP)
    return build_ellipsoidal_backward_certifier(
        problem,
        cfg;
        symbolic_builder = symbolic_builder,
    )
end

function main(cfg::DoublePendulumMPPIConfig = DoublePendulumMPPIConfig())
    Random.seed!(cfg.rng_seed)

    input_mapping = build_double_pendulum_input_mapping(cfg)

    generator_builder = (problem, system_cfg, control_cfg, cfg_) ->
        build_double_pendulum_mppi_generator(
            problem,
            system_cfg,
            control_cfg,
            cfg_;
            input_mapping = input_mapping,
        )

    certifier_builder = (problem, system_cfg, control_cfg, cfg_) ->
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
        build_concrete_system = () -> build_double_pendulum_system_cfg(cfg; pendulum_module = DP),
        build_control_problem = () -> build_double_pendulum_control_cfg(cfg; pendulum_module = DP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )

    warmup_candidate = OP.get_seed(run_result.solver.generator)
    warmup_candidate === nothing && error("MPPI warmup trajectory is missing.")

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
