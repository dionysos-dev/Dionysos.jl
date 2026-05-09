include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pendulum", "double_pendulum.jl"))
const DP = DoublePendulum

Base.@kwdef struct DoublePendulumBenchmarkConfig
    l1::Float64 = 1.0
    l2::Float64 = 1.0
    m1::Float64 = 1.0
    m2::Float64 = 1.0
    g::Float64 = 9.81
    objective::String = "benchmark_up_convex"

    Δt::Float64 = 0.1
    hx::SVector{4, Float64} = SVector(5 * (pi / 180), 5 * (pi / 180), 0.25, 0.25)
    periodic_dims::SVector{2, Int} = SVector(1, 2)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 300
    trajectory_mode::Symbol = :abstract_traj
    input_values::Tuple{Vararg{Float64}} = vals_tuple = Tuple(-3.5:0.25:3.5)

    terminal_radius::Float64 = 1.5
    λ::Float64 = 0.05
    maxδx::Float64 = 120.0
    maxδu::Float64 = 60.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
    )
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.1, 0.1), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    plot_mp4::Bool = true
    verbose::Bool = false
end

function build_double_pendulum_generator(
    problem,
    _system_cfg,
    control_cfg,
    cfg::DoublePendulumBenchmarkConfig;
    input_mapping,
)
    return build_centered_generator(
        problem,
        cfg;
        input_mapping = input_mapping,
        jacobian_bound = DP.jacobian_bound(;
            l1 = cfg.l1,
            l2 = cfg.l2,
            m1 = cfg.m1,
            m2 = cfg.m2,
            g = cfg.g,
            ωmax = 6.0,
        ),
        x0_provider = _ -> control_cfg.x0,
        trajectory_mode = cfg.trajectory_mode,
    )
end

function build_double_pendulum_certifier(
    problem,
    _system_cfg,
    _control_cfg,
    cfg::DoublePendulumBenchmarkConfig,
)
    symbolic_builder = build_double_pendulum_symbolic_builder(cfg; pendulum_module = DP)
    return build_ellipsoidal_backward_certifier(
        problem,
        cfg;
        symbolic_builder = symbolic_builder,
    )
end

function main(cfg::DoublePendulumBenchmarkConfig = DoublePendulumBenchmarkConfig())
    input_mapping = build_double_pendulum_input_mapping(cfg)

    generator_builder = (problem, system_cfg, control_cfg, cfg_) ->
        build_double_pendulum_generator(
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
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            angles_title = "double_pendulum - (theta1, theta2)",
            velocities_title = "double_pendulum - (omega1, omega2)",
            state_title = "double_pendulum - states",
            control_title = "double_pendulum - control",
            phase_title = "double_pendulum - phase portraits",
        )
        animation_paths = save_double_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            run_result.nominal_candidate;
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            l1 = cfg.l1,
            l2 = cfg.l2,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "double_pendulum",
        )
        return (; gif = animation_paths.gif_path, mp4 = animation_paths.mp4_path)
    end

    return run_benchmark(
        cfg;
        scenario_name = "double_pendulum",
        build_concrete_system = () -> build_double_pendulum_system_cfg(cfg; pendulum_module = DP),
        build_control_problem = () -> build_double_pendulum_control_cfg(cfg; pendulum_module = DP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
