# premier test du 2trailors
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

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
Configuration for the marche arriere benchmark.
"""
Base.@kwdef struct MarcheArriereConfig
    Δt::Float64 = 0.2
    hx::SVector{5, Float64} =
        SVector(0.6, 0.6, 6 * (pi / 180), 8 * (pi / 180), 8 * (pi / 180))
    periodic_dims::SVector{3, Int} = SVector(3, 4, 5)
    periodic_periods::SVector{3, Float64} = SVector(2pi, 2pi, 2pi)
    periodic_start::SVector{3, Float64} = SVector(-pi, -pi, -pi)
    nstep::Int = 300

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
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi, -pi),
        SVector(19.0, 10.0, pi, pi, pi),
    )
    x_domain =
        AV.with_phi_limits(x_domain; phi1_max = deg2rad(65.0), phi2_max = deg2rad(65.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(8.0, -1.0), SVector(19.0, 4.0)),
        UT.HyperRectangle(SVector(8.0, 7.0), SVector(19.0, 10.0)),
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
        [2.0, -0.5],
        [2.0, 0.5],
        [-2.0, -0.5],
        [-2.0, 0.5],
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
        SVector(17.5, 5.0, pi - deg2rad(6.0), -deg2rad(8.0), -deg2rad(8.0)),
        SVector(19.0, 6.4, pi + deg2rad(6.0), deg2rad(8.0), deg2rad(8.0)),
    )

    return (; x0, initial_set, target_set)
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    return run_vehicle_benchmark(
        cfg;
        scenario_name = "double_remorque_marche_avant",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
