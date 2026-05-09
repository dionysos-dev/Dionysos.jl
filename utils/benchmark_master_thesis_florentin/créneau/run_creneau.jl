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
Base.@kwdef struct MarcheArriereConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
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

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    return run_vehicle_benchmark(
        cfg;
        scenario_name = "creneau",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
    )
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
