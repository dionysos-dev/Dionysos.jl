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
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 40

    terminal_radius::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    symbolic_rk4_substeps::Int = 1
    ΔX::IA.IntervalBox{4, Float64} = IA.IntervalBox(
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.1, 0.1),
        IA.interval(-0.1, 0.1),
    )
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    verbose::Bool = false

    seed_trajectory_mode::Symbol = :abstract_traj
    seed_num_substeps::Int = 5
    mppi_nsamples::Int = 1900
    mppi_niter::Int = 5 * 15
    mppi_λ::Float64 = 0.3
    mppi_noise_v::Float64 = 0.15
    mppi_noise_σ::Float64 = 0.08
end


######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(0.0, 0.0, -pi, -pi),
        SVector(21.0,5.6, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(65.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(0.0, 0.0), SVector(6.3, 2.1)), # voiture arrière
        UT.HyperRectangle(SVector(14.3, 0.0), SVector(21.0, 2.1)), # voiture avant
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = 0.959931 
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 2.2, L2 = 2.8, Lc = 0.9) # on devra peut etre rendre cela plus réaliste
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping()
    inputs_delta = [
        [-2.0, 0.0], # vitesse sans angle
        [0.0, 0.0],
        [2.0, 0.0],
        [-1.0, 0.1],[-1.0, -0.1],
        [1.0, 0.1],[1.0, -0.1],
        

        [-1.0, 0.35],[-1.0, -0.35],
        [1.0, 0.35],[1.0, -0.35],

        [-1.0, 0.5],[-1.0, -0.5],
    ]

    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    # départ dans la voie, devant la voiture avant
    x0 = SVector(19.0, 3.8, 0.0, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(18.0, 3.5, -deg2rad(3.0), -deg2rad(2.0)),
        SVector(20.6, 4.5,  deg2rad(3.0),  deg2rad(2.0)),
    )

    target_set = UT.HyperRectangle(
        SVector(10.0, 0.5, -deg2rad(6.0), -deg2rad(3.0)),
        SVector(12.8, 1.6,  deg2rad(6.0),  deg2rad(3.0)),
    )

    return (; x0, initial_set, target_set)
end

function build_wrap_state(cfg::MarcheArriereConfig)
    return ST.get_periodic_wrapper(
        cfg.periodic_dims,
        cfg.periodic_periods;
        start = cfg.periodic_start,
    )
end

function periodic_state_error(x, xref, cfg::MarcheArriereConfig)
    e = collect(Float64, x .- xref)

    for i in eachindex(cfg.periodic_dims)
        d = Int(cfg.periodic_dims[i])
        p = Float64(cfg.periodic_periods[i])
        e[d] = mod(e[d] + 0.5 * p, p) - 0.5 * p
    end

    return SVector{4, Float64}(e)
end

function project_input_to_domain(u, u_domain)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function build_mppi_generator(
    problem,
    system_cfg,
    control_cfg,
    cfg::MarcheArriereConfig;
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
    target_center = UT.get_center(control_cfg.target_set)

    discrete_dynamics = (prob, x, u, k, Δt) -> f_disc(x, u)

    noise_sampler = function (rng, u, k)
        return [
            cfg.mppi_noise_v * Random.randn(rng),
            cfg.mppi_noise_σ * Random.randn(rng),
        ]
    end

    get_reference_states = function ()
        seed_cand = OP.get_trajectory(seed_gen)
        seed_cand === nothing && return nothing
        return ST.enum_elems(seed_cand.x_traj)
    end

    project_input = u -> project_input_to_domain(u, system_cfg.u_domain)

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
        scenario_name = "creneau_mppi",
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
        title12 = "creneau_mppi (x,y)",
        title34 = "creneau_mppi (theta,phi)",
    )
    stat_result = run_kappa_statistical_check(run_result; n_samples = 500)
    save_kappa_statistical_plots!(stat_result; wrap_angles = true)

    return run_result
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
