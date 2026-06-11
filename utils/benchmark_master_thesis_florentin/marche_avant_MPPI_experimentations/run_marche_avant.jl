# c'est fait à la va vite, il faut fixe ça demain dans la nouvelle version de dio
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
import LinearAlgebra as LA
import MathematicalSystems as MS
import Random

# Load the articulated vehicle benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the marche avant benchmark.
"""
Base.@kwdef struct MarcheAvantConfig
    Δt::Float64 = 0.1
    hx::SVector{4, Float64} = SVector(0.2, 0.2, 3 * (pi / 180), 3 * (pi / 180)) # ref baseline
    periodic_dims::SVector{1, Int} = SVector(3)
    periodic_periods::SVector{1, Float64} = SVector(2pi)
    periodic_start::SVector{1, Float64} = SVector(-pi)
    nstep::Int = 60

    terminal_radius::Float64 = 0.45
    use_terminal_john_ellipsoid::Bool = true
    terminal_john_shrink::Float64 = 1.0
    terminal_success_distance2::Float64 = 0.5
    terminal_outside_weight::Float64 = 1.0e4
    terminal_center_weight::Float64 = 100.0
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    state_scaling_mode::Symbol = :matrix # :none, :std, or :matrix
    state_scaling_std_floor::Float64 = 1.0e-3
    # In :matrix mode, only the diagonal is used as the coordinate scaling.
    state_scaling_matrix::Matrix{Float64} = Matrix{Float64}(LA.I, 4, 4)
    symbolic_rk4_substeps::Int = 1

    # Boites de linearisation adaptatives, exprimees en coordonnees physiques.
    # Elles sont ensuite utilisees de facon coherente avec le scaling interne.
    adaptive_linearization_boxes::Bool = true
    ΔX_initial::Vector{Float64} = [0.10, 0.10, 0.05, 0.05]
    ΔX_min::Vector{Float64} = [1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4]
    ΔX_max::Vector{Float64} = [12.0, 12.0, 2pi, 2deg2rad(55.0)]
    ΔU_initial::Vector{Float64} = [0.10, 0.10]
    ΔU_min::Vector{Float64} = [1.0e-4, 1.0e-4]
    ΔU_max::Vector{Float64} = [10.0, 3.0]
    adaptive_box_growth::Float64 = 1.5
    adaptive_box_safety::Float64 = 1.05
    adaptive_box_max_iters::Int = 30
    adaptive_box_atol::Float64 = 1.0e-8
    adaptive_box_verbose::Bool = false
    adaptive_box_search_scales::Vector{Float64} = [0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 3.0]
    adaptive_box_objective::Symbol = :max_volume
    adaptive_box_keep_first_consistent::Bool = false
    strict_reach_avoid_audit::Bool = true
    strict_reach_avoid_audit_verbose::Bool = true

    ΔX::IA.IntervalBox{4, Float64} =
        IA.IntervalBox(
            IA.interval(-0.2, 0.2), # on peut aller jusqu'à -1.0 ; 1.0
            IA.interval(-0.2, 0.2),
            IA.interval(-0.2, 0.2),
            IA.interval(-0.2, 0.2),
        ) * 5
    ΔU::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-0.2, 0.2), IA.interval(-0.2, 0.2)) * 5
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
    mppi_nsamples::Int = 1800
    mppi_niter::Int = 15
    mppi_λ::Float64 = 0.3

    # Echelle du bruit ajoute sur les deux composantes du controle.
    # Ici on garde quelque chose de simple et lisible : un bruit gaussien
    # independant sur la vitesse et sur la commande de braquage.
    mppi_noise_v::Float64 = 0.45
    mppi_noise_σ::Float64 = 0.15

    # Verification statistique post-certification : on sample dans
    # l'ellipsoide initial certifie et on rejoue la chaine de feedbacks.
    kappa_statistical_samples::Int = 0
    kappa_statistical_seed::Int = 1

    # Stress-test empirique post-certification.  Pour alpha > 1, les points
    # sont hors garantie formelle : on mesure seulement le bassin de succes en
    # simulation sous la chaine de feedback certifiee sur E0.
    run_empirical_inflation_stress_test::Bool = false
    empirical_inflation_factors::Vector{Float64} =
        [1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0]
    empirical_inflation_samples_per_alpha::Int = 5000
    empirical_inflation_shell_sampling::Bool = false
    empirical_inflation_seed::Int = 2
    empirical_inflation_rollout_alpha::Float64 = 2.0
    empirical_inflation_save_rollout_plot::Bool = false
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(SVector(-1.0, -1.0, -pi, -pi), SVector(10.0, 10.0, pi, pi))
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(75.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 7.0), SVector(10.0, 10.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = 0.959931 # ce qui est donné dans l'article
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 2.2, L2 = 2.8, Lc = 0.9)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system, obstacles_xy)
end

function distance_to_xy_obstacle(x, obstacle::UT.HyperRectangle)
    dx = max(obstacle.lb[1] - x[1], 0.0, x[1] - obstacle.ub[1])
    dy = max(obstacle.lb[2] - x[2], 0.0, x[2] - obstacle.ub[2])
    return hypot(dx, dy)
end

function min_distance_to_xy_obstacles(x, obstacles_xy)
    isempty(obstacles_xy) && return Inf
    return minimum(obstacle -> distance_to_xy_obstacle(x, obstacle), obstacles_xy)
end

function build_terminal_john_ellipsoid(lower, upper; shrink = 1.0)
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
        throw(ArgumentError("terminal_john_shrink must be finite and in (0, 1]"))

    radii = 0.5 .* (ub .- lb)
    center = 0.5 .* (lb .+ ub)
    shape = Matrix(LA.Diagonal(1.0 ./ (shrink .* radii) .^ 2))
    LA.isposdef(LA.Symmetric(shape)) ||
        throw(ArgumentError("terminal John ellipsoid shape must be positive definite"))

    return center, shape
end

function build_terminal_john_ellipsoid(target_set; shrink = 1.0)
    hasproperty(target_set, :lb) && hasproperty(target_set, :ub) ||
        throw(ArgumentError("target_set must expose lb and ub bounds"))
    return build_terminal_john_ellipsoid(target_set.lb, target_set.ub; shrink = shrink)
end

function terminal_john_ellipsoid_data(problem, cfg::MarcheAvantConfig)
    cfg.use_terminal_john_ellipsoid || return nothing
    terminal_center, terminal_shape =
        build_terminal_john_ellipsoid(problem.target_set; shrink = cfg.terminal_john_shrink)
    return (; terminal_center, terminal_shape)
end

function terminal_ellipsoidal_distance2(x, terminal_center, terminal_shape)
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

function trajectory_state_std_scaling(candidate; floor::Float64 = 1.0e-3)
    xs = collect(ST.enum_elems(candidate.x_traj))
    isempty(xs) && return fill(1.0, 4)

    nx = length(first(xs))
    n = length(xs)
    means = [sum(x[i] for x in xs) / n for i in 1:nx]

    if n == 1
        return fill(1.0, nx)
    end

    stds = [sqrt(sum((x[i] - means[i])^2 for x in xs) / (n - 1)) for i in 1:nx]
    return max.(stds, floor)
end

function matrix_state_scaling(M)
    # A = Matrix{Float64}(M)
    # size(A, 1) == size(A, 2) || error("state_scaling_matrix must be square")
    # return [abs(A[i, i]) for i in 1:size(A, 1)]
    return [10.0, 10.0, deg2rad(55.0), deg2rad(55.0)]*1.5
end # [2.0, 1.5, 0.35, 0.35] # [0.75, 0.60, 0.25, 0.25] #[10.0, 10.0, deg2rad(55.0), deg2rad(55.0)] 

function resolve_state_scaling(cfg::MarcheAvantConfig, candidate)
    cfg.state_scaling_mode == :none && return nothing
    cfg.state_scaling_mode == :std &&
        return trajectory_state_std_scaling(candidate; floor = cfg.state_scaling_std_floor)
    cfg.state_scaling_mode == :matrix &&
        return matrix_state_scaling(cfg.state_scaling_matrix)
    return error(
        "Unknown state_scaling_mode $(cfg.state_scaling_mode). Use :none, :std, or :matrix.",
    )
end

function build_input_mapping() # c'est pour la génération de la trajectoire
    inputs_delta = [
        [0.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
        [2.0, -0.55],
        [2.0, 0.55],
        [1.0, 0.0],
        [-1.0, 0.0],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    initial_set = UT.HyperRectangle(
        SVector(-0.4, -0.4, -20 * (pi / 180), -20 * (pi / 180)),
        SVector(0.4, 0.4, 20 * (pi / 180), 20 * (pi / 180)),
    )
    target_set = UT.HyperRectangle(
        SVector(8.5, 5.0, -35 * (pi / 180), -35 * (pi / 180)),
        SVector(10.0, 6.7, 35 * (pi / 180), 35 * (pi / 180)),
    )
    return (; x0, initial_set, target_set)
end

# MPPI travaille avec des trajectoires concretes.
# On reutilise donc le meme wrapping periodique que dans le reste de la
# pipeline pour que les etats restent dans la representation attendue.
function build_wrap_state(cfg::MarcheAvantConfig)
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
function periodic_state_error(x, xref, cfg::MarcheAvantConfig)
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
    cfg::MarcheAvantConfig;
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
    target_center = UT.get_center(problem.target_set)
    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
    obstacles_xy = system_cfg.obstacles_xy

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    # Bruit additif simple, volontairement sans sophistication.
    noise_sampler = function (rng, u, k)
        return [cfg.mppi_noise_v * Random.randn(rng), cfg.mppi_noise_σ * Random.randn(rng)]
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)
    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return ST.enum_elems(seed_cand.x_traj)
    end

    # trajectory_cost = function (prob, cand)
    #     xs = collect(ST.enum_elems(cand.x_traj))
    #     us = collect(ST.enum_elems(cand.u_traj))
    #     return
    # end 

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

        # lissage des commandes
        w_u = 0.03
        w_du = 0.25
        w_ddu = 0.05

        # coût de marge / proximité frontière
        w_margin = 400000.0
        clearance_radius = 1.2

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

            terminal_hit = if terminal_data === nothing
                xw ∈ prob.target_set
            else
                terminal_ellipsoidal_distance2(
                    xw,
                    terminal_data.terminal_center,
                    terminal_data.terminal_shape,
                ) <= cfg.terminal_success_distance2 + 1.0e-8
            end

            if terminal_hit
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
            clearance = min_distance_to_xy_obstacles(xw, obstacles_xy)
            if clearance < clearance_radius
                J += w_margin * ((clearance_radius - clearance) / clearance_radius)^2
            end
        end

        # -----------------------------------------------------
        # 5) coût terminal
        # -----------------------------------------------------
        xT = wrap_state(xs[min(hit_index, length(xs))])
        eT = periodic_state_error(xT, target_center, cfg)

        J += w_goal_pos_T * (eT[1]^2 + eT[2]^2)
        J += w_goal_th_T * (eT[3]^2)
        J += w_goal_phi_T * (eT[4]^2)

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
        (prob, x) -> wrap_state(x);
        truncate_at_first_target = terminal_data === nothing,
        postprocess_candidate = postprocess_candidate,
    )

    return OP.MPPIGenerator(problem, mppi_cfg, seed_gen)
end

function _build_marche_avant_mppi_certifier(
    problem,
    system_cfg,
    cfg::MarcheAvantConfig;
    vehicle_module,
    symbolic_builder = build_symbolic_builder(vehicle_module, system_cfg.params),
    backend = build_backend(; verbose = cfg.verbose),
    state_scaling = nothing,
)
    terminal_data = terminal_john_ellipsoid_data(problem, cfg)
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
        state_scaling = state_scaling,
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
        backend,
        opts,
    )
    return SC.EllipsoidalBackwardCertifier(cert_cfg, symbolic_builder)
end

mutable struct MarcheAvantScalingCertifier <: SC.AbstractSymbolicCertifier
    problem::Any
    system_cfg::Any
    cfg::MarcheAvantConfig
    vehicle_module::Any
    symbolic_builder::Any
    backend::Any
    candidate::Any
    inner::Any
    result::Any
    success::Bool
    solve_time_sec::Float64
end

function build_marche_avant_mppi_certifier(
    problem,
    system_cfg,
    cfg::MarcheAvantConfig;
    vehicle_module,
    symbolic_builder = build_symbolic_builder(vehicle_module, system_cfg.params),
    backend = build_backend(; verbose = cfg.verbose),
)
    return MarcheAvantScalingCertifier(
        problem,
        system_cfg,
        cfg,
        vehicle_module,
        symbolic_builder,
        backend,
        nothing,
        nothing,
        nothing,
        false,
        0.0,
    )
end

function SC.set_problem!(
    cert::MarcheAvantScalingCertifier,
    problem::Dionysos.Problem.ProblemType,
)
    cert.problem = problem
    cert.candidate = nothing
    cert.inner = nothing
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function SC.set_trajectory!(cert::MarcheAvantScalingCertifier, candidate)
    cert.candidate = candidate
    cert.inner = nothing
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function SC.certify!(cert::MarcheAvantScalingCertifier)
    state_scaling = resolve_state_scaling(cert.cfg, cert.candidate)
    cert.inner = _build_marche_avant_mppi_certifier(
        cert.problem,
        cert.system_cfg,
        cert.cfg;
        vehicle_module = cert.vehicle_module,
        symbolic_builder = cert.symbolic_builder,
        backend = cert.backend,
        state_scaling = state_scaling,
    )
    SC.set_problem!(cert.inner, cert.problem)
    SC.set_trajectory!(cert.inner, cert.candidate)
    SC.certify!(cert.inner)

    cert.result = SC.get_result(cert.inner)
    cert.success = SC.get_success(cert.inner)
    cert.solve_time_sec = SC.get_solve_time(cert.inner)

    if cert.success && cert.cfg.strict_reach_avoid_audit
        audit = audit_reach_avoid_ellipsoids(cert.problem, cert.system_cfg, cert.result)
        if !audit.success
            cert.cfg.strict_reach_avoid_audit_verbose && println(
                "reach-avoid audit failed: ",
                length(audit.violations),
                " ellipsoid/domain-obstacle violations; first = ",
                isempty(audit.violations) ? nothing : first(audit.violations),
            )
            cert.success = false
            cert.result =
                certification_result_with_reach_avoid_failure(cert.result, audit.failed_k)
        end
    end

    return cert
end

SC.get_result(cert::MarcheAvantScalingCertifier) = cert.result
SC.get_success(cert::MarcheAvantScalingCertifier) = cert.success
SC.get_solve_time(cert::MarcheAvantScalingCertifier) = cert.solve_time_sec

function _safe_solve_time(obj)
    obj === nothing && return NaN
    try
        return Float64(OP.get_solve_time(obj))
    catch
        try
            return Float64(SC.get_solve_time(obj))
        catch
            return NaN
        end
    end
end

function marche_avant_timing_summary(run_result)
    gen = run_result.solver.generator
    seed_gen = hasproperty(gen, :seed_generator) ? gen.seed_generator : nothing
    cert = run_result.solver.certifier

    abstract_trajectory_sec = _safe_solve_time(seed_gen)
    mppi_with_seed_sec = _safe_solve_time(gen)
    mppi_sec =
        isfinite(mppi_with_seed_sec) && isfinite(abstract_trajectory_sec) ?
        max(0.0, mppi_with_seed_sec - abstract_trajectory_sec) : mppi_with_seed_sec
    lmi_sec = _safe_solve_time(cert)
    total_sec = _safe_solve_time(run_result.solver)
    total_components_sec =
        sum(t for t in (abstract_trajectory_sec, mppi_sec, lmi_sec) if isfinite(t))

    return (;
        abstract_trajectory_sec,
        mppi_sec,
        mppi_with_seed_sec,
        lmi_sec,
        total_sec,
        total_components_sec,
    )
end

function print_marche_avant_timing_summary(timings)
    println("\n=== Marche avant timing summary ===")
    println("abstract_trajectory_sec = ", timings.abstract_trajectory_sec)
    println("mppi_sec                = ", timings.mppi_sec)
    println("mppi_with_seed_sec      = ", timings.mppi_with_seed_sec)
    println("lmi_sec                 = ", timings.lmi_sec)
    println("total_sec               = ", timings.total_sec)
    println("total_components_sec    = ", timings.total_components_sec)
    return nothing
end

function save_marche_avant_timing_summary(path::AbstractString, timings)
    open(path, "w") do io
        println(io, "metric,seconds")
        println(io, "abstract_trajectory_sec,", timings.abstract_trajectory_sec)
        println(io, "mppi_sec,", timings.mppi_sec)
        println(io, "mppi_with_seed_sec,", timings.mppi_with_seed_sec)
        println(io, "lmi_sec,", timings.lmi_sec)
        println(io, "total_sec,", timings.total_sec)
        return println(io, "total_components_sec,", timings.total_components_sec)
    end
    return path
end

function active_state_scaling(cert::MarcheAvantScalingCertifier)
    cert.inner === nothing && return nothing
    return cert.inner.config.options.state_scaling
end

function ellipsoid_axis_aligned_bounds(E)
    c = collect(Float64, E.c)
    Q = LA.inv(LA.Symmetric(Matrix{Float64}(E.P)))
    radii = sqrt.(max.(LA.diag(Q), 0.0))
    return c .- radii, c .+ radii
end

function state_domain_base_set(X)
    return X isa UT.LazySetMinus ? X.A : X
end

function ellipsoid_inside_hyperrectangle(E, rect; tol::Float64 = 1.0e-8)
    lo, hi = ellipsoid_axis_aligned_bounds(E)
    n = min(length(lo), length(rect.lb))
    lower_violation = maximum(Float64(rect.lb[i]) - lo[i] for i in 1:n)
    upper_violation = maximum(hi[i] - Float64(rect.ub[i]) for i in 1:n)
    violation = max(lower_violation, upper_violation)
    return violation <= tol, violation
end

function ellipsoid_xy_aabb_intersects_obstacle(E, obstacle; tol::Float64 = 1.0e-8)
    lo, hi = ellipsoid_axis_aligned_bounds(E)
    return !(
        hi[1] < Float64(obstacle.lb[1]) - tol ||
        lo[1] > Float64(obstacle.ub[1]) + tol ||
        hi[2] < Float64(obstacle.lb[2]) - tol ||
        lo[2] > Float64(obstacle.ub[2]) + tol
    )
end

function certified_ellipsoid_state_index(i::Int, nsteps::Int)
    # lmi_data.ellipsoids is stored as [E_terminal, E_K, E_{K-1}, ..., E_1].
    return nsteps + 2 - i
end

function audit_reach_avoid_ellipsoids(
    problem,
    system_cfg,
    cert_result;
    tol::Float64 = 1.0e-8,
)
    cert_result.success || return (;
        success = false,
        failed_k = cert_result.failed_k,
        violations = NamedTuple[],
    )
    cert_result.lmi_data === nothing &&
        return (; success = false, failed_k = nothing, violations = NamedTuple[])
    !hasproperty(cert_result.lmi_data, :ellipsoids) &&
        return (; success = false, failed_k = nothing, violations = NamedTuple[])

    Xbase = state_domain_base_set(problem.system.X)
    obstacles_xy = hasproperty(system_cfg, :obstacles_xy) ? system_cfg.obstacles_xy : ()
    ellipsoids = cert_result.lmi_data.ellipsoids
    nsteps = length(cert_result.steps)
    violations = NamedTuple[]

    for (i, E) in enumerate(ellipsoids)
        k_state = certified_ellipsoid_state_index(i, nsteps)
        inside_domain, domain_violation =
            ellipsoid_inside_hyperrectangle(E, Xbase; tol = tol)
        if !inside_domain
            push!(
                violations,
                (;
                    kind = :domain_containment,
                    ellipsoid_index = i,
                    state_index = k_state,
                    violation = domain_violation,
                ),
            )
        end

        for (j, obstacle) in enumerate(obstacles_xy)
            if ellipsoid_xy_aabb_intersects_obstacle(E, obstacle; tol = tol)
                push!(
                    violations,
                    (;
                        kind = :obstacle_aabb_intersection,
                        ellipsoid_index = i,
                        state_index = k_state,
                        obstacle_index = j,
                        violation = Inf,
                    ),
                )
            end
        end
    end

    failed_k = isempty(violations) ? nothing : min(first(violations).state_index, nsteps)
    return (; success = isempty(violations), failed_k, violations)
end

function certification_result_with_reach_avoid_failure(cert_result, failed_k)
    return SC.EllipsoidalCertificationResult(
        false,
        failed_k,
        cert_result.solve_time_sec,
        cert_result.steps,
        cert_result.controller,
        cert_result.lmi_data,
    )
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

function main(cfg::MarcheAvantConfig = MarcheAvantConfig())
    run_result = run_vehicle_benchmark(
        cfg;
        scenario_name = "marche_avant_mppi",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
        generator_builder = build_mppi_generator,
        certifier_builder = build_marche_avant_mppi_certifier,
    )
    timings = marche_avant_timing_summary(run_result)
    print_marche_avant_timing_summary(timings)
    timing_csv = save_marche_avant_timing_summary(
        joinpath(run_result.outputs.root, "timings.csv"),
        timings,
    )
    run_result = merge(run_result, (; timings, timing_csv))

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
        title12 = "marche_avant_mppi (x,y)",
        title34 = "marche_avant_mppi (theta,phi)",
    )

    stat_result = nothing
    if cert_result.success && cfg.kappa_statistical_samples > 0
        stat_result = run_kappa_statistical_check(
            run_result;
            n_samples = cfg.kappa_statistical_samples,
            num_substeps = cfg.seed_num_substeps,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            project_inputs = false,
            verbose = true,
            rng = Random.MersenneTwister(cfg.kappa_statistical_seed),
        )
        save_kappa_statistical_plots!(
            stat_result;
            output_dir = plots_dir,
            basename = "kappa_statistical_rollouts",
            wrap_angles = true,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            show_ellipsoids = false,
            show_exit_points = true,
            axis_labels_12 = (L"x_1\,[\mathrm{m}]", L"x_2\,[\mathrm{m}]"),
            axis_labels_34 = (L"\theta_1\,[\mathrm{rad}]", L"\phi\,[\mathrm{rad}]"),
        )
    elseif cfg.kappa_statistical_samples > 0
        println("Skipping kappa statistical check because certification failed.")
    end

    empirical_stress_result = nothing
    if cert_result.success && cfg.run_empirical_inflation_stress_test
        empirical_stress_result = run_empirical_inflation_stress_test(
            run_result;
            inflation_factors = cfg.empirical_inflation_factors,
            samples_per_alpha = cfg.empirical_inflation_samples_per_alpha,
            shell_sampling = cfg.empirical_inflation_shell_sampling,
            selected_rollout_alpha = cfg.empirical_inflation_rollout_alpha,
            num_substeps = cfg.seed_num_substeps,
            periodic_dims = cfg.periodic_dims,
            periodic_periods = cfg.periodic_periods,
            periodic_start = cfg.periodic_start,
            project_inputs = false,
            verbose = true,
            rng = Random.MersenneTwister(cfg.empirical_inflation_seed),
        )
        stress_csv = save_empirical_inflation_stress_csv(
            joinpath(run_result.outputs.root, "empirical_inflation_stress.csv"),
            empirical_stress_result.summary_rows,
        )
        stress_plot = plot_empirical_inflation_stress_rates(
            empirical_stress_result.summary_rows;
            output_dir = plots_dir,
            filename = "empirical_success_vs_inflation.pdf",
        )

        stress_rollout_paths = nothing
        if cfg.empirical_inflation_save_rollout_plot &&
           empirical_stress_result.selected_rollout_result !== nothing
            alpha_tag = replace(string(cfg.empirical_inflation_rollout_alpha), "." => "_")
            stress_rollout_paths = save_kappa_statistical_plots!(
                empirical_stress_result.selected_rollout_result;
                output_dir = plots_dir,
                basename = "empirical_inflation_rollouts_alpha_$(alpha_tag)",
                wrap_angles = true,
                periodic_dims = cfg.periodic_dims,
                periodic_periods = cfg.periodic_periods,
                periodic_start = cfg.periodic_start,
                show_ellipsoids = false,
                show_exit_points = true,
                axis_labels_12 = (L"x_1\,[\mathrm{m}]", L"x_2\,[\mathrm{m}]"),
                axis_labels_34 = (L"\theta_1\,[\mathrm{rad}]", L"\phi\,[\mathrm{rad}]"),
                rollout_lw = 0.45,
                rollout_alpha = 0.18,
                initial_sample_alpha = 0.34,
            )
        end

        empirical_stress_result = merge(
            empirical_stress_result,
            (;
                outputs = (;
                    csv = stress_csv,
                    plot = stress_plot,
                    rollout_paths = stress_rollout_paths,
                )
            ),
        )
    elseif cfg.run_empirical_inflation_stress_test
        println("Skipping empirical inflation stress test because certification failed.")
    end

    return merge(
        run_result,
        (; statistical = stat_result, empirical_inflation_stress = empirical_stress_result),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
