include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "pendulum",
        "simple_pendulum.jl",
    ),
)
const SP = SimplePendulum

Base.@kwdef struct SimplePendulumBenchmarkConfig
    l::Float64 = 1.0
    g::Float64 = 9.81
    objective::String = "benchmark_up_convex"

    Δt::Float64 = 0.05
    hx::SVector{2, Float64} = SVector(5 * (pi / 180), 0.25)
    periodic_dims::SVector{1, Int} = SVector(1)
    periodic_periods::SVector{1, Float64} = SVector(2pi)
    periodic_start::SVector{1, Float64} = SVector(-pi)
    nstep::Int = 6600
    trajectory_mode::Symbol = :abstract_traj
    input_values::Tuple{Vararg{Float64}} = Tuple(-2.3:0.25:2.3)

    terminal_radius::Float64 = 0.35
    λ::Float64 = 0.05
    maxδx::Float64 = 80.0
    maxδu::Float64 = 40.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-0.05, 0.05), IA.interval(-0.05, 0.05))
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.05, 0.05), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    plot_mp4::Bool = true
    verbose::Bool = false
end

function build_simple_pendulum_generator(
    problem,
    _system_cfg,
    control_cfg,
    cfg::SimplePendulumBenchmarkConfig;
    input_mapping,
)
    return build_centered_generator(
        problem,
        cfg;
        input_mapping = input_mapping,
        jacobian_bound = SP.jacobian_bound(cfg.l, cfg.g),
        x0_provider = _ -> control_cfg.x0,
        trajectory_mode = cfg.trajectory_mode,
    )
end

function build_simple_pendulum_certifier(
    problem,
    _system_cfg,
    _control_cfg,
    cfg::SimplePendulumBenchmarkConfig,
)
    symbolic_builder = build_pendulum_symbolic_builder(cfg; pendulum_module = SP)
    return build_ellipsoidal_backward_certifier(
        problem,
        cfg;
        symbolic_builder = symbolic_builder,
    )
end

function main(cfg::SimplePendulumBenchmarkConfig = SimplePendulumBenchmarkConfig())
    input_mapping = build_pendulum_input_mapping(cfg)

    generator_builder =
        (problem, system_cfg, control_cfg, cfg_) -> build_simple_pendulum_generator(
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
        save_pendulum_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            run_result.nominal_candidate;
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            phase_title = "simple_pendulum - phase space",
            state_title = "simple_pendulum - states",
            control_title = "simple_pendulum - control",
        )
        animation_paths = save_pendulum_animation!(
            run_result.outputs.animations_dir,
            run_result.problem,
            run_result.nominal_candidate;
            cert_result = run_result.result.certification,
            show_ellipsoids = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            save_gif = cfg.plot_gif,
            save_mp4 = cfg.plot_mp4,
            title = "simple_pendulum",
        )
        return (; gif = animation_paths.gif_path, mp4 = animation_paths.mp4_path)
    end

    return run_benchmark(
        cfg;
        scenario_name = "simple_pendulum",
        build_concrete_system = () -> build_pendulum_system_cfg(cfg; pendulum_module = SP),
        build_control_problem = () -> build_pendulum_control_cfg(cfg; pendulum_module = SP),
        generator_builder = generator_builder,
        certifier_builder = certifier_builder,
        save_artifacts! = save_artifacts!,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
