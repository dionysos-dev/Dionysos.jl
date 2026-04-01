# c'est encore un benchmark "script" a l'ancienne :
# on reutilise la pipeline existante telle quelle, en l'instanciant
# proprement pour le nouveau probleme flat-plate glider 7D.
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import MathematicalSystems as MS
import Random

# Load the flat-plate glider benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "flat_plat_glider_7D.jl"))
const FG = FlatPlateGlider7D

######################################################
###################[DATA & INIT]######################
######################################################

Base.@kwdef struct FlatPlateGliderConfig
    params::FG.Params{Float64} = FG.Params()

    Δt::Float64 = 0.05
    hx::SVector{7, Float64} = SVector(
        0.50,                # px
        0.50,                # pz
        6 * (pi / 180),     # θ
        6 * (pi / 180),     # φ
        2.0,                 # vx
        2.0,                 # vz
        1.0,                 # q
    )
    periodic_dims::SVector{0, Int} = SVector{0, Int}()
    periodic_periods::SVector{0, Float64} = SVector{0, Float64}()
    periodic_start::SVector{0, Float64} = SVector{0, Float64}()
    nstep::Int = 30

    terminal_radius::Float64 = 0.30
    λ::Float64 = 0.001
    maxδx::Float64 = 150.0
    maxδu::Float64 = 150.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{7, Float64} = IA.IntervalBox(IA.interval(-0.01, 0.01), 7)
    ΔU::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(-0.01, 0.01), 1)
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    # V1 "survie" : on coupe le GIF pour garder un benchmark leger.
    plot_gif::Bool = false
    verbose::Bool = true

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 2
    mppi_nsamples::Int = 120
    mppi_niter::Int = 8
    mppi_λ::Float64 = 2.0
    mppi_noise_phidot::Float64 = 1.5

    # Le papier laisse q(tf) non borne ; ici on fixe une fenetre finie
    # pour rester compatible avec les HyperRectangles et le certifier.
    target_q_abs_bound::Float64 = 6.0

    perch_pos::SVector{2, Float64} = SVector(0.0, 0.0)
    perch_radius::Float64 = 0.05
    plot_margin_x::Float64 = 0.25
    plot_margin_z::Float64 = 0.25
end

function pipeline_only_cfg(cfg::FlatPlateGliderConfig)
    return FlatPlateGliderConfig(
        params = cfg.params,
        Δt = cfg.Δt,
        hx = cfg.hx,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        nstep = cfg.nstep,
        terminal_radius = cfg.terminal_radius,
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        output_root = cfg.output_root,
        plot_subdir = cfg.plot_subdir,
        animation_subdir = cfg.animation_subdir,
        plot_gif = false,
        verbose = cfg.verbose,
        seed_trajectory_mode = cfg.seed_trajectory_mode,
        seed_num_substeps = cfg.seed_num_substeps,
        mppi_nsamples = cfg.mppi_nsamples,
        mppi_niter = cfg.mppi_niter,
        mppi_λ = cfg.mppi_λ,
        mppi_noise_phidot = cfg.mppi_noise_phidot,
        target_q_abs_bound = cfg.target_q_abs_bound,
        perch_pos = cfg.perch_pos,
        perch_radius = cfg.perch_radius,
        plot_margin_x = cfg.plot_margin_x,
        plot_margin_z = cfg.plot_margin_z,
    )
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_sets(cfg::FlatPlateGliderConfig)
    p = cfg.params

    # V1 de benchmark simplifie :
    # on reste sur un probleme local de reachability, mais avec un
    # domaine plus large et une cible moins triviale que le simple
    # "vol tout droit" de la V0.
    x_domain = UT.HyperRectangle(
        SVector(
            1.5,
            -1.0,
            -pi / 2,
            p.phi_min,
            0.0,
            -4.0,
            -8.0,
        ),
        SVector(
            6.0,
            1.0,
            pi / 2,
            p.phi_max,
            10.0,
            2.0,
            8.0,
        ),
    )
    x_domain = FG.with_phi_limit(x_domain; phi_min = p.phi_min, phi_max = p.phi_max)

    u_domain = UT.HyperRectangle(SVector(-p.phi_dot_max), SVector(p.phi_dot_max))

    # Etat nominal de l'article.
    x0 = SVector(3.5, 0.1, 0.0, 0.0, 7.0, 0.0, 0.0)

    # On remet une petite incertitude initiale non negligeable pour
    # garder un vrai probleme abstrait.
    initial_set = UT.HyperRectangle(
        x0 - SVector(0.5, 0.05, 0.05, 0.02, 0.2, 0.2, 0.5),
        x0 + SVector(0.5, 0.05, 0.05, 0.02, 0.2, 0.2, 0.5),
    )

    # Cible toujours large et alignee avec la grille coarse, mais plus
    # proche de l'esprit "approche / orientation / ralentissement".
    qTf = min(cfg.target_q_abs_bound, p.q_abs_bound)
    target_set = UT.HyperRectangle(
        SVector(
            4.5,
            -0.5,
            10 * (pi / 180),
            p.phi_min,
            2.0,
            -2.5,
            -qTf,
        ),
        SVector(
            5.5,
            0.5,
            70 * (pi / 180),
            p.phi_max,
            8.0,
            0.5,
            qTf,
        ),
    )

    return (; x_domain, u_domain, x0, initial_set, target_set)
end

function build_concrete_system(cfg::FlatPlateGliderConfig)
    sets = build_sets(cfg)
    concrete_system = FG.system(sets.x_domain; _U_ = sets.u_domain, params = cfg.params)
    return (; x_domain = sets.x_domain, u_domain = sets.u_domain, params = cfg.params, concrete_system)
end

function build_control_problem(cfg::FlatPlateGliderConfig)
    sets = build_sets(cfg)
    return (; x0 = sets.x0, initial_set = sets.initial_set, target_set = sets.target_set)
end

function build_input_mapping(cfg::FlatPlateGliderConfig)
    p = cfg.params
    inputs = [
        [-p.phi_dot_max],
        [-0.5 * p.phi_dot_max],
        [0.0],
        [0.5 * p.phi_dot_max],
        [p.phi_dot_max],
    ]
    return MP.ListMapping(inputs)
end

function abstraction_periodicity(cfg::FlatPlateGliderConfig)
    isempty(cfg.periodic_dims) && return nothing
    return periodicity_kwargs(cfg)
end

function build_wrap_state(cfg::FlatPlateGliderConfig)
    isempty(cfg.periodic_dims) && return identity
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

function periodic_state_error(x, xref, cfg::FlatPlateGliderConfig)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    return SVector{7, Float64}(e)
end

function project_input_to_domain(u, u_domain)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function build_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::FlatPlateGliderConfig;
    vehicle_module,
    input_mapping,
)
    seed_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        input_mapping,
        vehicle_module.jacobian_bound(system_cfg.params),
        abstraction_periodicity(cfg),
        cfg.nstep,
        _ -> control_cfg.x0;
        num_substeps = cfg.seed_num_substeps,
        trajectory_mode = cfg.seed_trajectory_mode,
    )

    seed_gen = OP.CenteredAbstractionGenerator(problem, seed_cfg)

    disc_system = ST.discretize_continuous_system(
        system_cfg.concrete_system,
        cfg.Δt;
        num_substeps = cfg.seed_num_substeps,
    )
    f_disc = MS.mapping(disc_system)
    wrap_state = build_wrap_state(cfg)
    target_center = UT.get_center(control_cfg.target_set)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    noise_sampler = function (rng, u, k)
        return [cfg.mppi_noise_phidot * Random.randn(rng)]
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)

    trajectory_cost = function (prob, cand)
        xs = collect(ST.enum_elems(cand.x_traj))
        us = collect(ST.enum_elems(cand.u_traj))

        BAD_COST = 1.0e15
        J = 0.0

        # Poids simples de debug / benchmark.
        w_step = 3.0
        w_pos = 12.0
        w_theta = 2.0
        w_phi = 0.4
        w_vel = 0.8
        w_q = 0.3
        w_u = 0.01
        w_du = 0.08
        w_ddu = 0.02

        hit_target = false
        hit_index = length(xs)

        for k in eachindex(xs)
            xw = wrap_state(xs[k])

            if !(xw ∈ prob.system.X)
                return BAD_COST
            end

            if xw ∈ prob.target_set
                hit_target = true
                hit_index = k
                J /= 20.0
                break
            end

            e = periodic_state_error(xw, target_center, cfg)
            J += w_step
            J += w_pos * (e[1]^2 + e[2]^2)
            J += w_theta * e[3]^2
            J += w_phi * e[4]^2
            J += w_vel * (e[5]^2 + e[6]^2)
            J += w_q * e[7]^2
        end

        if !hit_target
            xT = wrap_state(last(xs))
            eT = periodic_state_error(xT, target_center, cfg)
            J += 800.0 * (eT[1]^2 + eT[2]^2)
            J += 120.0 * eT[3]^2
            J += 40.0 * eT[4]^2
            J += 50.0 * (eT[5]^2 + eT[6]^2)
            J += 25.0 * eT[7]^2
            J += 1.0e4
        end

        last_u_index = min(length(us), max(hit_index - 1, 0))
        for k in 1:last_u_index
            u = us[k]
            J += w_u * sum(abs2, u)

            if k >= 2
                du = u .- us[k - 1]
                J += w_du * sum(abs2, du)
            end

            if k >= 3
                ddu = u .- 2 .* us[k - 1] .+ us[k - 2]
                J += w_ddu * sum(abs2, ddu)
            end
        end

        return J
    end

    success_fun = function (prob, cand)
        xT = wrap_state(last(ST.enum_elems(cand.x_traj)))
        return xT ∈ prob.target_set
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
        (prob, x) -> wrap_state(x),
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function save_named_state_space_plots!(
    plots_dir::AbstractString,
    basename::AbstractString,
    problem,
    candidate;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims,
    periodic_periods,
    periodic_start,
    title12::AbstractString,
    title34::AbstractString,
)
    mktempdir() do tmp_dir
        save_state_space_plots!(
            tmp_dir,
            problem,
            candidate;
            cert_result = cert_result,
            show_ellipsoids = show_ellipsoids,
            unwrap_angles = unwrap_angles,
            wrap_angles = wrap_angles,
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
            periodic_start = periodic_start,
            title12 = title12,
            title34 = title34,
        )

        cp(
            joinpath(tmp_dir, "state_space_12.pdf"),
            joinpath(plots_dir, "$(basename)_12.pdf");
            force = true,
        )
        cp(
            joinpath(tmp_dir, "state_space_34.pdf"),
            joinpath(plots_dir, "$(basename)_34.pdf");
            force = true,
        )
    end

    return nothing
end

function save_glider_rollout!(
    vehicle_module,
    concrete_system,
    params,
    cfg::FlatPlateGliderConfig,
    x_traj,
    u_traj;
    giffile,
    fps = 20,
    every = 1,
    dt = cfg.Δt,
    title = nothing,
)
    gr()

    xl = (
        concrete_system.X.lb[1] - cfg.plot_margin_x,
        concrete_system.X.ub[1] + cfg.plot_margin_x,
    )
    zl = (
        concrete_system.X.lb[2] - cfg.plot_margin_z,
        concrete_system.X.ub[2] + cfg.plot_margin_z,
    )

    draw_params = vehicle_module.DrawParams(params)

    return vehicle_module.live_glider_progression(
        params,
        draw_params,
        x_traj,
        u_traj,
        xl,
        zl;
        domain = concrete_system.X,
        every = every,
        dt = dt,
        giffile = giffile,
        fps = fps,
        title = title,
        perch_pos = cfg.perch_pos,
        perch_radius = cfg.perch_radius,
    )
end

function main(cfg::FlatPlateGliderConfig = FlatPlateGliderConfig())
    run_result = run_vehicle_benchmark(
        pipeline_only_cfg(cfg);
        scenario_name = "flat_plat_glider_mppi",
        vehicle_module = FG,
        build_concrete_system = () -> build_concrete_system(cfg),
        build_control_problem = () -> build_control_problem(cfg),
        input_mapping = build_input_mapping(cfg),
        generator_builder = build_mppi_generator,
    )
    warmup_candidate = OP.get_seed(run_result.solver.generator)
    warmup_candidate === nothing && error("MPPI warmup trajectory is missing.")

    plots_dir = run_result.outputs.plots_dir
    problem = run_result.problem
    cert_result = run_result.result.certification

    # On reecrit explicitement les plots principaux avec des titres x-z,
    # car le helper historique a ete pense pour des coordonnees x-y.
    save_named_state_space_plots!(
        plots_dir,
        "state_space",
        problem,
        run_result.nominal_candidate;
        cert_result = cert_result,
        show_ellipsoids = true,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "flat_plat_glider_mppi (x,z)",
        title34 = "flat_plat_glider_mppi (theta,phi)",
    )

    save_named_state_space_plots!(
        plots_dir,
        "abstract_traj_used_as_warmup",
        problem,
        warmup_candidate;
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "Abstract traj used as warmup",
        title34 = "Abstract traj used as warmup",
    )

    save_named_state_space_plots!(
        plots_dir,
        "mppi_candidate_traj",
        problem,
        run_result.nominal_candidate;
        cert_result = nothing,
        show_ellipsoids = false,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "MPPI candidate traj (x,z)",
        title34 = "MPPI candidate traj (theta,phi)",
    )

    save_named_state_space_plots!(
        plots_dir,
        "ellipsoids",
        problem,
        run_result.nominal_candidate;
        cert_result = cert_result,
        show_ellipsoids = true,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "flat_plat_glider_mppi certified tube (x,z)",
        title34 = "flat_plat_glider_mppi certified tube (theta,phi)",
    )

    gif_path = nothing
    if cfg.plot_gif
        gif_path = joinpath(run_result.outputs.animations_dir, "rollout.gif")
        save_glider_rollout!(
            FG,
            problem.system,
            cfg.params,
            cfg,
            run_result.nominal_candidate.x_traj,
            run_result.nominal_candidate.u_traj;
            giffile = gif_path,
            dt = run_result.nominal_candidate.Ts,
            title = "flat_plat_glider_mppi - pipeline certifie",
        )
        println("gif = ", gif_path)
    end

    return merge(run_result, (; gif = gif_path))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
