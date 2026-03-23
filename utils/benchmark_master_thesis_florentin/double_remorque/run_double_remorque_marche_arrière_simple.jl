# premier test du 2trailors
include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

# Load the articulated vehicle benchmark model from Dionysos problems.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle_2trailers.jl"))
const AV = ArticulatedVehicle2Trailers

######################################################
###################[DATA & INIT]######################
######################################################
"""
Configuration for the marche arriere benchmark.
"""
Base.@kwdef struct MarcheArriereConfig
    Δt::Float64 = 0.2
    hx::SVector{5, Float64} = SVector(0.25, 0.25, 2.5 * (pi / 180), 2.5 * (pi / 180), 2.5 * (pi / 180))
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
        IA.interval(-0.05, 0.05),
        IA.interval(-0.05, 0.05),
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
        SVector(0.0, 3.0, -pi, -pi, -pi),
        SVector(8.0, 5.0,  pi, pi,  pi),
    )
    x_domain = AV.with_phi_limits(
        x_domain;
        phi1_max = deg2rad(50.0),
        phi2_max = deg2rad(50.0),
    )

    obstacles_xy = Any[]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    # marche arrière sans braquage
    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 1.0, L2 = 1.0, L3 = 1.0, Lc = 0.5, Lc2 = 0.0)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)

    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping() # c'est pour la génération de la trajectoire
    inputs_delta = [
        [-1.0, 0.0],
        [0.0, 0.0],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end 

function build_control_problem()
    x0 = SVector(7.5, 4.0, 0.0, 0.0, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(7.0, 3.45, -deg2rad(6.0), -deg2rad(6.0), -deg2rad(6.0)),
        SVector(8.0, 4.55,  deg2rad(6.0),  deg2rad(6.0),  deg2rad(6.0)),
    )

    target_set = UT.HyperRectangle(
        SVector(0.0, 3.0, -deg2rad(8.0), -deg2rad(8.0), -deg2rad(8.0)),
        SVector(1.0, 5.0,  deg2rad(8.0),  deg2rad(8.0),  deg2rad(8.0)),
    )

    return (; x0, initial_set, target_set)
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    return run_vehicle_benchmark(
        cfg;
        scenario_name = "double_remorque_marche_arriere",
        vehicle_module = AV,
        build_concrete_system = build_concrete_system,
        build_control_problem = build_control_problem,
        input_mapping = build_input_mapping(),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
