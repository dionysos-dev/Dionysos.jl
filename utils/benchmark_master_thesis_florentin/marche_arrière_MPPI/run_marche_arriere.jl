# c'est fait à la va vite, il faut fixe ça demain dans la nouvelle version de dio
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import MathematicalSystems as MS
import Random

# Load the articulated vehicle benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the marche arriere benchmark.
"""
Base.@kwdef struct MarcheArriereConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 3 * (pi / 180), 2 * (pi / 180))
    #hx::SVector{4, Float64} = SVector(0.5, 0.5, 5 * (pi / 180), 5 * (pi / 180)) # en closed loop je peux faire une discrétisation largement moins fine et trouver une traj malgès tout! (bon pour MPPI)
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 60

    terminal_radius::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.2, 0.2),
        IA.interval(-0.2, 0.2),
        IA.interval(-0.2, 0.2),
        IA.interval(-0.2, 0.2),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.2, 0.2),
        IA.interval(-0.2, 0.2),
    )
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    verbose::Bool = false

    # Parametres de la seed fournie a MPPI.
    # On conserve un generateur centred abstraction comme source de
    # warm-start, mais MPPI re-simule ensuite toujours sa trajectoire
    # propre avec la dynamique concrete.
    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 5

    # Parametres MPPI separes de ceux du certifier.
    # On ne reutilise pas `cfg.λ` car ce champ existe deja dans la
    # pipeline de certification et porte un autre sens.
    mppi_nsamples::Int = 50000
    mppi_niter::Int = 4*15
    mppi_λ::Float64 = 0.3

    # Echelle du bruit ajoute sur les deux composantes du controle.
    # Ici on garde quelque chose de simple et lisible : un bruit gaussien
    # independant sur la vitesse et sur la commande de braquage.
    mppi_noise_v::Float64 = 0.35
    mppi_noise_σ::Float64 = 0.15
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi),
        SVector(16.0, 10.0, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(55.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(9.0, -1.0), SVector(16.0, 4.0)),
        UT.HyperRectangle(SVector(9.0, 7.4), SVector(16.0, 10.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = 0.959931 # ce qui est donné dans l'article
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 2.2, L2 = 2.8, Lc = 0.9)
    
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping() # c'est pour la génération de la trajectoire
    inputs_delta = [
        [2.0, -0.5],
        [2.0, 0.5],
        [-2.0, -0.5],
        [-2.0, 0.5],

        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],

        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    initial_set = UT.HyperRectangle(
        SVector(-1.0, -1.0, -0.4, -0.4),
        SVector(1.0, 1.0, 0.4, 0.4),
    )
    target_set = UT.HyperRectangle(
        SVector(15.0, 5.0, pi - 3 * (pi / 180), -2 * (pi / 180)),
        SVector(16.0, 6.2, pi + 3 * (pi / 180),  2 * (pi / 180)),
    )
    return (; x0, initial_set, target_set)
end

# MPPI travaille avec des trajectoires concretes.
# On reutilise donc le meme wrapping periodique que dans le reste de la
# pipeline pour que les etats restent dans la representation attendue.
function build_wrap_state(cfg::MarcheArriereConfig)
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

# Pour comparer un etat au centre de la cible, on calcule l'erreur sur
# les dimensions periodiques en choisissant toujours le plus petit ecart
# modulo la periode. Sans cela, un angle proche de -pi et un angle proche
# de pi sembleraient artificiellement tres loin.
function periodic_state_error(x, xref, cfg::MarcheArriereConfig)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    return SVector{4, Float64}(e)
end

# MPPI manipule librement des controles perturbes.
# On reprojette donc explicitement chaque controle dans le domaine
# admissible du systeme avant de le simuler.
function project_input_to_domain(u, u_domain)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

# Constructeur local du generateur MPPI pour ce benchmark.
# On garde toute la pipeline commune inchangée ; seule la facon de
# construire le generateur nominal varie.
function build_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::MarcheArriereConfig;
    vehicle_module,
    input_mapping,
)
    # La seed est produite par le generateur historique centered
    # abstraction. On s'en sert uniquement pour recuperer une premiere
    # sequence de controles `u_traj`.
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

    #on devrait réutiliser celle def dans symbolic problem ici 
    disc_system = ST.discretize_continuous_system(
        system_cfg.concrete_system,
        cfg.Δt;
        num_substeps = cfg.seed_num_substeps,
    )
    f_disc = MS.mapping(disc_system)
    wrap_state = build_wrap_state(cfg)
    target_center = UT.get_center(control_cfg.target_set)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    # Bruit additif simple, volontairement sans sophistication.
    noise_sampler = function (rng, u, k)
        return [
            cfg.mppi_noise_v * Random.randn(rng),
            cfg.mppi_noise_σ * Random.randn(rng),
        ]
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)
    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return ST.enum_elems(seed_cand.x_traj)
    end
    # Cout volontairement simple :
    # - rapprocher la trajectoire de la cible,
    # - penaliser moderement l'effort de commande,
    # - penaliser fortement les sorties du domaine admissible.
    #
    # Le but de cette V1 n'est pas d'encoder une fonction de cout
    # parfaite, mais d'avoir quelque chose de clair et directement
    # exploitable pour lancer le benchmark.
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
        w_step        = 1.0

        # suivi d'une trajectoire de référence (ou nominale)
        w_ref_pos     = 2.0
        w_ref_ang     = 0.25

        # coût terminal vers la cible
        w_goal_pos_T  = 800.0
        w_goal_th_T   = 120.0
        w_goal_phi_T  = 120.0
        w_miss        = 1.0e4

        # lissage des commandes # je devrais modifier ici afin de rendre mes commandes moins wiggly
        w_u           = 0.03*3
        w_du          = 0.25*2
        w_ddu         = 0.05*3

        # coût de marge / proximité frontière
        w_margin      = 40.0

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
        J += w_goal_th_T  * (eT[3]^2)
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

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    run_result = run_vehicle_benchmark(
        cfg;
        scenario_name = "marche_arriere_mppi",
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
        title12 = "marche_arriere_mppi (x,y)",
        title34 = "marche_arriere_mppi (theta,phi)",
    )
    stat_result = run_kappa_statistical_check(run_result; n_samples = 500)
    save_kappa_statistical_plots!(stat_result; wrap_angles = true)
    return run_result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
