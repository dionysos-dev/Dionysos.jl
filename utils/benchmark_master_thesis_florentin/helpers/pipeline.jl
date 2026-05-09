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

function resolve_problem(system_cfg, control_cfg)
    if control_cfg !== nothing && hasproperty(control_cfg, :problem)
        prob = control_cfg.problem
        prob !== nothing && return prob
    end
    return build_problem(system_cfg, control_cfg)
end

function default_prepare_for_certification(cfg)
    if hasproperty(cfg, :periodic_dims) && hasproperty(cfg, :periodic_periods)
        return build_periodic_certification_preprocessor(
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
        )
    end
    return identity
end

function build_centered_generator(
    problem,
    cfg;
    Δt = cfg.Δt,
    hx = cfg.hx,
    input_mapping,
    jacobian_bound,
    x0_provider,
    periodicity = periodicity_kwargs(cfg),
    nstep = cfg.nstep,
    num_substeps::Int = 5,
    trajectory_mode::Symbol = :abstract_traj,
)
    gen_cfg = OP.CenteredAbstractionConfig(
        Δt,
        hx,
        input_mapping,
        jacobian_bound,
        periodicity,
        nstep,
        x0_provider;
        num_substeps = num_substeps,
        trajectory_mode = trajectory_mode,
    )

    return OP.CenteredAbstractionGenerator(problem, gen_cfg)
end

function build_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg;
    vehicle_module,
    input_mapping,
)
    return build_centered_generator(
        problem,
        cfg;
        input_mapping,
        jacobian_bound = vehicle_module.jacobian_bound(system_cfg.params),
        x0_provider = _ -> control_cfg.x0,
        trajectory_mode = :abstract_traj, # :closed_loop # :abstract_traj
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

function build_ellipsoidal_backward_certifier(
    problem,
    cfg;
    symbolic_builder,
    backend = build_backend(; verbose = cfg.verbose),
    terminal_center = nothing,
    terminal_shape = nothing,
)
    state_scaling = if hasproperty(cfg, :use_state_scaling) && !cfg.use_state_scaling
        nothing
    elseif hasproperty(cfg, :state_scaling)
        raw = cfg.state_scaling
        raw === nothing ? nothing : Float64.(collect(raw))
    else
        nothing
    end

    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        terminal_center = terminal_center,
        terminal_shape = terminal_shape,
        state_scaling = state_scaling,
    )

    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        backend,
        opts,
    )

    return SC.EllipsoidalBackwardCertifier(cert_cfg, symbolic_builder)
end

function build_certifier(
    problem,
    system_cfg,
    cfg;
    vehicle_module,
    symbolic_builder = build_symbolic_builder(vehicle_module, system_cfg.params),
    backend = build_backend(; verbose = cfg.verbose),
)
    return build_ellipsoidal_backward_certifier(
        problem,
        cfg;
        symbolic_builder = symbolic_builder,
        backend = backend,
    )
end

function _report_pipeline_status(gen, cert, solver, nominal_candidate, cert_candidate, paths, gif_path, mp4_path)
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
    mp4_path !== nothing && println("mp4 = ", mp4_path)
    return nothing
end

function run_benchmark(
    cfg;
    scenario_name::AbstractString,
    build_concrete_system,
    build_control_problem,
    generator_builder,
    certifier_builder,
    save_artifacts! = nothing,
    prepare_for_certification = default_prepare_for_certification(cfg),
)
    paths = output_paths(cfg)

    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = resolve_problem(system_cfg, control_cfg)

    gen = generator_builder(problem, system_cfg, control_cfg, cfg)
    cert = certifier_builder(problem, system_cfg, control_cfg, cfg)
    solver = OP.CertifiedPipelineSolver(gen, cert)

    OP.set_problem!(solver, problem)
    OP.solve!(solver; prepare_for_certification = prepare_for_certification)

    result = OP.get_result(solver)
    result === nothing && error("Pipeline returned no result.")
    result.candidate === nothing && error("Pipeline produced no candidate trajectory.")
    result.certification_candidate === nothing &&
        error("Pipeline produced no certification trajectory.")

    nominal_candidate = result.candidate
    cert_candidate = result.certification_candidate

    run_result = (
        ;
        solver,
        result,
        problem,
        config = cfg,
        outputs = paths,
        system_cfg,
        control_cfg,
        nominal_candidate,
        certification_candidate = cert_candidate,
    )

    extras = if save_artifacts! === nothing
        (;)
    else
        saved = save_artifacts!(run_result)
        saved === nothing ? (;) : saved
    end

    gif_path = hasproperty(extras, :gif) ? extras.gif : nothing
    mp4_path = hasproperty(extras, :mp4) ? extras.mp4 : nothing
    _report_pipeline_status(
        gen,
        cert,
        solver,
        nominal_candidate,
        cert_candidate,
        paths,
        gif_path,
        mp4_path,
    )

    return merge(run_result, extras)
end

function run_vehicle_benchmark(
    cfg;
    scenario_name::AbstractString,
    vehicle_module,
    build_concrete_system,
    build_control_problem,
    input_mapping,
    # On laisse injectable la construction du generateur nominal.
    # Par defaut, on garde le comportement historique du benchmark
    # via `build_generator`. Mais ce petit point d'extension permet
    # d'utiliser un autre heuristic generator, par exemple MPPI,
    # sans dupliquer toute la pipeline ni modifier les helpers aval.
    generator_builder = build_generator,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
)
    wrapped_generator_builder = function (problem, system_cfg, control_cfg, cfg)
        return generator_builder(
            problem,
            system_cfg,
            control_cfg,
            cfg;
            vehicle_module = vehicle_module,
            input_mapping = input_mapping,
        )
    end

    wrapped_certifier_builder = function (problem, system_cfg, _control_cfg, cfg)
        return build_certifier(
            problem,
            system_cfg,
            cfg;
            vehicle_module = vehicle_module,
        )
    end

    save_artifacts! = function (run_result)
        save_state_space_plots!(
            run_result.outputs.plots_dir,
            run_result.problem,
            run_result.nominal_candidate;
            cert_result = run_result.result.certification,
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
            gif_path = joinpath(run_result.outputs.animations_dir, "rollout.gif")
            plot_articulated_vehicle!(
                vehicle_module,
                run_result.problem.system,
                run_result.system_cfg.params,
                run_result.nominal_candidate.x_traj,
                run_result.nominal_candidate.u_traj;
                giffile = gif_path,
                dt = run_result.nominal_candidate.Ts,
                title = "$(scenario_name) - pipeline certifie",
            )
        end

        return (; gif = gif_path)
    end

    return run_benchmark(
        cfg;
        scenario_name = scenario_name,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        generator_builder = wrapped_generator_builder,
        certifier_builder = wrapped_certifier_builder,
        save_artifacts! = save_artifacts!,
    )
end
