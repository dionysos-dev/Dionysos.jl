import MathOptInterface as MOI
import StaticArrays: SVector
using JuMP: optimizer_with_attributes
using MosekTools

periodicity_kwargs(cfg) = (;
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

function output_paths(cfg)
    root = cfg.output_root
    plots_dir = joinpath(root, cfg.plot_subdir)
    animations_dir = joinpath(root, cfg.animation_subdir)
    mkpath(plots_dir)
    mkpath(animations_dir)
    return (; root, plots_dir, animations_dir)
end

function build_backend(; verbose::Bool = false)
    return optimizer_with_attributes(
        MosekTools.Optimizer,
        MOI.Silent() => !verbose,
        MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => (verbose ? 10 : 0),
        MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
    )
end

function build_problem(system_cfg, control_cfg)
    return PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

function build_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg;
    vehicle_module,
    input_mapping,
)
    gen_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        input_mapping,
        vehicle_module.jacobian_bound(system_cfg.params),
        periodicity_kwargs(cfg),
        cfg.nstep,
        _ -> control_cfg.x0,
    )

    return OP.CenteredAbstractionGenerator{
        typeof(problem),
        typeof(gen_cfg),
        Any,
        Any,
    }(
        problem,
        gen_cfg,
        nothing,
        nothing,
        false,
        0.0,
    )
end

function build_symbolic_builder(vehicle_module, params)
    return function (prob, candidate, certifier_cfg)
        o = certifier_cfg.options
        return vehicle_module.symbolic_system(
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

function build_certifier(
    problem,
    system_cfg,
    cfg;
    vehicle_module,
    symbolic_builder = build_symbolic_builder(vehicle_module, system_cfg.params),
    backend = build_backend(; verbose = cfg.verbose),
)
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

    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        backend,
        opts,
    )

    return SC.EllipsoidalBackwardCertifier{
        typeof(problem),
        Any,
        typeof(cert_cfg),
        Any,
        typeof(symbolic_builder),
    }(
        nothing,
        nothing,
        cert_cfg,
        nothing,
        false,
        0.0,
        symbolic_builder,
    )
end

function _report_pipeline_status(gen, cert, solver, nominal_candidate, cert_candidate, paths, gif_path)
    println("generator_success = ", OP.get_success(gen))
    println("certifier_success = ", SC.get_success(cert))
    println("pipeline_success = ", OP.get_success(solver))
    println("candidate_horizon = ", OP.horizon(nominal_candidate))
    println("certification_candidate_horizon = ", OP.horizon(cert_candidate))

    result = OP.get_result(solver)
    if result !== nothing && result.certification !== nothing && hasproperty(result.certification, :steps)
        println("cert_steps = ", length(result.certification.steps))
        println("failed_k = ", result.certification.failed_k)
    end

    println("output_root = ", paths.root)
    println("plots_dir = ", paths.plots_dir)
    gif_path !== nothing && println("gif = ", gif_path)
    return nothing
end

function run_vehicle_benchmark(
    cfg;
    scenario_name::AbstractString,
    vehicle_module,
    build_concrete_system,
    build_control_problem,
    input_mapping,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
)
    paths = output_paths(cfg)

    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)

    gen = build_generator(
        problem,
        system_cfg,
        control_cfg,
        cfg;
        vehicle_module = vehicle_module,
        input_mapping = input_mapping,
    )
    cert = build_certifier(
        problem,
        system_cfg,
        cfg;
        vehicle_module = vehicle_module,
    )
    solver = OP.CertifiedPipelineSolver(gen, cert)

    prepare_for_certification = build_periodic_certification_preprocessor(
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
    )

    OP.set_problem!(solver, problem)
    OP.solve!(solver; prepare_for_certification = prepare_for_certification)

    result = OP.get_result(solver)
    result === nothing && error("Pipeline returned no result.")
    result.candidate === nothing && error("Pipeline produced no candidate trajectory.")
    result.certification_candidate === nothing &&
        error("Pipeline produced no certification trajectory.")

    nominal_candidate = result.candidate
    cert_candidate = result.certification_candidate

    save_state_space_plots!(
        paths.plots_dir,
        problem,
        nominal_candidate;
        cert_result = result.certification,
        show_ellipsoids = show_ellipsoids,
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "$(scenario_name) (x,y)",
        title34 = "$(scenario_name) (theta,phi)",
    )

    gif_path = nothing
    if cfg.plot_gif
        gif_path = joinpath(paths.animations_dir, "rollout.gif")
        plot_articulated_vehicle!(
            vehicle_module,
            problem.system,
            system_cfg.params,
            nominal_candidate.x_traj,
            nominal_candidate.u_traj;
            giffile = gif_path,
            dt = nominal_candidate.Ts,
            title = "$(scenario_name) - pipeline certifie",
        )
    end

    _report_pipeline_status(gen, cert, solver, nominal_candidate, cert_candidate, paths, gif_path)

    return (;
        solver,
        result,
        problem,
        config = cfg,
        outputs = paths,
        nominal_candidate,
        certification_candidate = cert_candidate,
        gif = gif_path,
    )
end
