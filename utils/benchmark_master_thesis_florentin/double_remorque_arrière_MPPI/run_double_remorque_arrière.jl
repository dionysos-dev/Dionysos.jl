# premier test du 2trailors
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import MathematicalSystems as MS
import Random

# Load the articulated vehicle benchmark model from Dionysos problems.
include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "articulated_vehicle_2trailers.jl",
    ),
)
const AV = ArticulatedVehicle2Trailers

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the double remorque marche avant simple benchmark.
"""
Base.@kwdef struct DoubleRemorqueMarcheAvantSimpleConfig
    Δt::Float64 = 0.2
    hx::SVector{5, Float64} =
        SVector(0.6, 0.6, 6 * (pi / 180), 8 * (pi / 180), 8 * (pi / 180))
    periodic_dims::SVector{3, Int} = SVector(3, 4, 5)
    periodic_periods::SVector{3, Float64} = SVector(2pi, 2pi, 2pi)
    periodic_start::SVector{3, Float64} = SVector(-pi, -pi, -pi)
    nstep::Int = 90 # 54 pour l'abstract

    terminal_radius::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{5, Float64} = IA.IntervalBox(
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
    )
    ΔU::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-0.1, 0.1), IA.interval(-0.1, 0.1))
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    verbose::Bool = false

    # Parametres de la seed fournie a MPPI.
    # On conserve le generateur historique centered abstraction comme
    # warm-start, puis MPPI resimule sa propre trajectoire concrete.
    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 5

    # Parametres MPPI separes de ceux du certifier.
    mppi_nsamples::Int = 150000
    mppi_niter::Int = 20
    mppi_λ::Float64 = 4.3
    mppi_noise_v::Float64 = 0.85
    mppi_noise_σ::Float64 = 0.58
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-4.0, -1.0, -pi, -pi, -pi),
        SVector(19.0, 15.0, pi, pi, pi),
    )
    x_domain =
        AV.with_phi_limits(x_domain; phi1_max = deg2rad(65.0), phi2_max = deg2rad(65.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(8.0, -1.0), SVector(19.0, 2.0)),
        UT.HyperRectangle(SVector(8.0, 11.7), SVector(19.0, 15.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    # marche arrière sans braquage
    δ_max = 0.959931
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 2.2, L2 = 2.8, L3 = 2.2, Lc = 0.9, Lc2 = 0.9)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping() # c'est pour la génération de la trajectoire
    inputs_delta = [
        [0.0, 0.0],
        [1.0, 0.0],
        [-1.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
        [2.0, -0.45],
        [2.0, 0.45],
        [-2.0, -0.45],
        [-2.0, 0.45],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    x0 = SVector(0.0, 0.0, 0.0, 0.0, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(-1.0, -1.0, -deg2rad(6.0), -deg2rad(8.0), -deg2rad(8.0)),
        SVector(0.2, 0.2, deg2rad(6.0), deg2rad(8.0), deg2rad(8.0)),
    )

    target_set = UT.HyperRectangle(
        SVector(16.5, 5.0, pi - deg2rad(6.0), -deg2rad(16.0), -deg2rad(16.0)),
        SVector(19.0, 9.4, pi + deg2rad(6.0), deg2rad(16.0), deg2rad(16.0)),
    )

    return (; x0, initial_set, target_set)
end

function build_wrap_state(cfg::DoubleRemorqueMarcheAvantSimpleConfig)
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

function periodic_state_error(x, xref, cfg::DoubleRemorqueMarcheAvantSimpleConfig)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    return SVector{5, Float64}(e)
end

function project_input_to_domain(u, u_domain)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function build_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::DoubleRemorqueMarcheAvantSimpleConfig;
    vehicle_module,
    input_mapping,
)
    seed_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        input_mapping,
        vehicle_module.jacobian_bound(system_cfg.params),
        periodicity_kwargs(cfg),
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
    target_center = UT.get_center(problem.target_set)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    noise_sampler = function (rng, u, k)
        return [cfg.mppi_noise_v * Random.randn(rng), cfg.mppi_noise_σ * Random.randn(rng)]
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)
    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return ST.enum_elems(seed_cand.x_traj)
    end

    trajectory_cost = function (prob, cand)
        xs = ST.enum_elems(cand.x_traj)
        us = ST.enum_elems(cand.u_traj)
        reference_states = get_reference_states()

        BAD_COST = 1.0e12
        reference_states === nothing && return BAD_COST

        J = 0.0

        # -----------------------------
        # poids principaux
        # -----------------------------
        w_step = 1.0

        # suivi d'une trajectoire de référence (ou nominale)
        w_ref_pos = 2.0
        w_ref_ang = 0.25

        # coût terminal vers la cible
        w_goal_pos_T = 800.0
        w_goal_th_T = 120.0
        w_goal_phi_T = 120.0
        w_miss = 1.0e4

        # lissage des commandes # je devrais modifier ici afin de rendre mes commandes moins wiggly
        w_u = 0.03*3
        w_du = 0.25*2
        w_ddu = 0.05*3

        # coût de marge / proximité frontière
        w_margin = 40.0

        # seuils de "handoff" façon Nav2
        near_goal_pos_radius = 1.5
        near_goal_ang_radius = 0.75

        hit_target = false
        hit_index = length(xs)

        for k in eachindex(xs)
            xw = wrap_state(xs[k])

            # -------------------------------------------------
            # 1) critic de faisabilité dure
            # -------------------------------------------------
            if !(xw ∈ prob.system.X)
                return BAD_COST
            end

            if xw ∈ prob.target_set
                hit_target = true
                hit_index = k
                break
            end

            # -------------------------------------------------
            # 2) critic de progression / suivi de référence
            # -------------------------------------------------
            # ref_x doit être la trajectoire nominale/abstraite alignée en temps
            ref_x = wrap_state(reference_states[min(k, length(reference_states))])

            e_ref = periodic_state_error(xw, ref_x, cfg)
            J += w_step
            J += w_ref_pos * (e_ref[1]^2 + e_ref[2]^2)
            J += w_ref_ang * (e_ref[3]^2 + e_ref[4]^2)

            # -------------------------------------------------
            # 3) critic goal activé seulement près du but
            # -------------------------------------------------
            e_goal = periodic_state_error(xw, target_center, cfg)
            dpos_goal = sqrt(e_goal[1]^2 + e_goal[2]^2)

            if dpos_goal <= near_goal_pos_radius
                J += 30.0 * (e_goal[1]^2 + e_goal[2]^2)
            end

            if dpos_goal <= near_goal_ang_radius
                J += 10.0 * e_goal[3]^2
                J += 10.0 * e_goal[4]^2
            end

            # -------------------------------------------------
            # 4) critic de marge de sécurité
            # -------------------------------------------------
            # Ce critic reste desactive tant qu'on ne dispose pas
            # d'une primitive de distance/coquille exploitable ici.
        end

        # -----------------------------------------------------
        # 5) coût terminal
        # -----------------------------------------------------
        xT = wrap_state(xs[min(hit_index, length(xs))])
        eT = periodic_state_error(xT, target_center, cfg)

        J += w_goal_pos_T * (eT[1]^2 + eT[2]^2)
        J += w_goal_th_T * (eT[3]^2)
        J += w_goal_phi_T * (eT[4]^2)

        if !hit_target
            J += w_miss
        end

        # -----------------------------------------------------
        # 6) critic sur les commandes
        # -----------------------------------------------------
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
        return any(x -> (wrap_state(x) ∈ prob.target_set), ST.enum_elems(cand.x_traj))
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
        return cp(
            joinpath(tmp_dir, "state_space_34.pdf"),
            joinpath(plots_dir, "$(basename)_34.pdf");
            force = true,
        )
    end

    return nothing
end

function main(
    cfg::DoubleRemorqueMarcheAvantSimpleConfig = DoubleRemorqueMarcheAvantSimpleConfig(),
)
    run_result = run_vehicle_benchmark(
        cfg;
        scenario_name = "double_remorque_marche_avant_mppi",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
        generator_builder = build_mppi_generator,
    )

    warmup_candidate = OP.get_seed(run_result.solver.generator)
    warmup_candidate === nothing && error("MPPI warmup trajectory is missing.")

    plots_dir = run_result.outputs.plots_dir
    problem = run_result.problem
    cert_result = run_result.result.certification

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
        title12 = "MPPI candidate traj",
        title34 = "MPPI candidate traj",
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
        title12 = "double_remorque_marche_avant_mppi (x,y)",
        title34 = "double_remorque_marche_avant_mppi (theta,phi1)",
    )
    stat_result = run_kappa_statistical_check(run_result; n_samples = 500)
    save_kappa_statistical_plots!(stat_result; wrap_angles = true)
    return run_result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
