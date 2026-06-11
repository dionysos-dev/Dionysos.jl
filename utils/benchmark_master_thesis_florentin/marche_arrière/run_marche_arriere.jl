# c'est fait à la va vite, il faut fixe ça demain dans la nouvelle version de dio
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

# Load the articulated vehicle benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the marche arriere benchmark.
"""
Base.@kwdef struct MarcheArriereConfig #dt = 0.1, nstep = 500
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    #hx::SVector{4, Float64} = SVector(0.5, 0.5, 5 * (pi / 180), 5 * (pi / 180)) # en abstract traj je peux faire une discrétisation largement moins fine et trouver une traj malgès tout! (bon pour MPPI)
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 300

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
    ΔU::IA.IntervalBox{2, Float64} =
        IA.IntervalBox(IA.interval(-0.2, 0.2), IA.interval(-0.2, 0.2))
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = true
    verbose::Bool = false
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(SVector(-1.0, -1.0, -pi, -pi), SVector(16.0, 10.0, pi, pi))
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(55.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(10.0, -1.0), SVector(16.0, 4.0)),
        UT.HyperRectangle(SVector(10.0, 8.0), SVector(16.0, 10.0)),
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
        # [2.0, -0.55],
        # [2.0, 0.55],
        # [-2.0, -0.55],
        # [-2.0, 0.55],

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
    initial_set =
        UT.HyperRectangle(SVector(-1.0, -1.0, -0.4, -0.4), SVector(1.0, 1.0, 0.4, 0.4))
    target_set = UT.HyperRectangle(
        SVector(13.5, 5.0, pi - 5 * (pi / 180), -5 * (pi / 180)),
        SVector(16.0, 6.2, pi + 5 * (pi / 180), 5 * (pi / 180)),
    )
    return (; x0, initial_set, target_set)
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    run_result = run_vehicle_benchmark(
        cfg;
        scenario_name = "marche_arriere",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
    )

    # je dois modifier la manière dont on sample et compte les succès
    # stat_result = run_kappa_statistical_check(
    #     run_result;
    #     n_samples = 300,
    #     num_substeps = 5,
    #     project_inputs = true,
    #     check_domain = true,
    # )

    # save_kappa_statistical_plots!(
    #     stat_result;
    #     wrap_angles = true,
    #     basename = "kappa_stats",
    # )

    # println(stat_result.summary)
    # println(stat_result.exit_histogram)
    return (; run_result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
