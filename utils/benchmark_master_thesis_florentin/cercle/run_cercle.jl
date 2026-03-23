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
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
        IA.interval(-0.15, 0.15),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.15, 0.15),
        IA.interval(-0.2, 0.2),
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
        SVector(12.0, 10.0,  pi,  pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(35.0))

    # obstacles: couloir en anneau (rectangle) autour d’un îlot central
    obstacles_xy = [
        UT.HyperRectangle(SVector(3.0, 4.0), SVector(9.0, 6.0)),  # mur gauche
        UT.HyperRectangle(SVector(4.0, 3.0), SVector(8.0, 7.0)),  # mur gauche
        UT.HyperRectangle(SVector(5.0, 2.0), SVector(7.0, 8.0)),  # mur gauche
    ]

    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system) 
end

function build_input_mapping() # c'est pour la génération de la trajectoire
    inputs_delta = [
        [-2.0, 0.0], # vitesse sans angle
        [0.0, 0.0],
        [-1.0, 0.0],
        [-0.5, 0.0],

        [-1.5, 0.1],[-1.5, -0.1], # moyen angles
        [-1.0, 0.1],[-1.0, -0.1],
        [-0.5, 0.1],[-0.5, -0.1],

        [-1.5, 0.35],[-1.5, -0.35],  # grand angles
        [-1.0, 0.35],[-1.0, -0.35],
        [-0.5, 0.35],[-0.5, -0.35],

        [-1.5, 0.65],[-1.5, -0.65], # très angles
        [-1.0, 0.65],[-1.0, -0.65],
        [-0.5, 0.65],[-0.5, -0.65],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    # départ côté est du couloir, orientation tangentielle
    x0 = SVector(6.0, 9.0, pi, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(5.25, 8.50, pi - deg2rad(4.0), -deg2rad(2.0)),
        SVector(6.45, 9.50, pi + deg2rad(4.0),  deg2rad(2.0)),
    )

    # même zone XY mais après un tour complet (theta décrémenté de 2π)
    target_set = UT.HyperRectangle(
        SVector(5.00, 0.50, - deg2rad(6.0), -deg2rad(4.0)),
        SVector(6.60, 2.0,  deg2rad(6.0),  deg2rad(4.0)),
    )

    return (; x0, initial_set, target_set)
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    return run_vehicle_benchmark(
        cfg;
        scenario_name = "cercle",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
