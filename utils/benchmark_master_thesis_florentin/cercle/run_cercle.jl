# c'est fait à la va vite, il faut fixe ça demain dans la nouvelle version de dio
import Dionysos
import IntervalArithmetic as IA
import MathOptInterface as MOI
import StaticArrays: SVector
using MosekTools
using JuMP

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const SY = DI.Symbolic
const SC = OP.SymbolicCertifier


######################################################
####################[INCLUDE]#########################
######################################################

include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
include(joinpath(@__DIR__, "..", "helpers", "plot.jl"))

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
        IA.interval(-0.5, 0.5),
        IA.interval(-0.5, 0.5),
        IA.interval(-0.5, 0.5),
        IA.interval(-0.5, 0.5),
    )
    ΔU::IA.IntervalBox{2, Float64} = IA.IntervalBox(
        IA.interval(-0.5, 0.5),
        IA.interval(-0.2, 0.2),
    )
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    output_root::String = joinpath(@__DIR__, "outputs")
    plot_subdir::String = "plots"
    animation_subdir::String = "animations"

    plot_gif::Bool = false
    verbose::Bool = false
end

periodicity_kwargs(cfg::MarcheArriereConfig) = (; # je n'aime pas l'idée qu'une structure soit stock pour wrap unwrap
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

function output_paths(cfg::MarcheArriereConfig)
    root = cfg.output_root
    plots_dir = joinpath(root, cfg.plot_subdir)
    animations_dir = joinpath(root, cfg.animation_subdir)
    mkpath(plots_dir)
    mkpath(animations_dir)
    return (; root, plots_dir, animations_dir)
end

function build_backend(; verbose::Bool = false)
    return optimizer_with_attributes(
        MosekTools.Optimizer,
        MOI.Silent() => !verbose,
        MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => (verbose ? 10 : 0),
        MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
    )
end

######################################################
############[INITIALISATION DU BENCHMARK]#############
######################################################

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi),
        SVector(10.0, 9.0, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(50.0))

    # obstacles: couloir en anneau (rectangle) autour d’un îlot central
    obstacles_xy = [
        UT.HyperRectangle(SVector(4.2, 3.2), SVector(5.8, 4.8)),  # ilot central
        UT.HyperRectangle(SVector(-1.0, -1.0), SVector(2.8, 9.0)), # mur gauche
        UT.HyperRectangle(SVector(7.2, -1.0), SVector(10.0, 9.0)), # mur droit
        UT.HyperRectangle(SVector(2.8, -1.0), SVector(7.2, 1.8)),  # mur bas
        UT.HyperRectangle(SVector(2.8, 6.2), SVector(7.2, 9.0)),   # mur haut
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
        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
        [1.0, -0.15],
        [1.0, 0.15],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return MP.ListMapping(inputs)
end

function build_control_problem()
    # départ côté est du couloir, orientation tangentielle
    x0 = SVector(6.6, 4.0, pi/2, 0.0)

    initial_set = UT.HyperRectangle(
        SVector(6.5, 3.9, pi/2 - 8*(pi/180), -5*(pi/180)),
        SVector(6.7, 4.1, pi/2 + 8*(pi/180),  5*(pi/180)),
    )

    # même zone XY mais après un tour complet (theta décrémenté de 2π)
    target_set = UT.HyperRectangle(
        SVector(6.5, 3.9, pi/2 - 2pi - 8*(pi/180), -5*(pi/180)),
        SVector(6.7, 4.1, pi/2 - 2pi + 8*(pi/180),  5*(pi/180)),
    )

    return (; x0, initial_set, target_set)
end

function build_problem(system_cfg, control_cfg)
    return PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

######################################################
##################[BUILD CERTIFIER]###################
######################################################

function build_generator(problem, system_cfg, control_cfg, cfg::MarcheArriereConfig)
    gen_cfg = OP.CenteredAbstractionConfig(
        cfg.Δt,
        cfg.hx,
        build_input_mapping(),
        AV.jacobian_bound(system_cfg.params),
        periodicity_kwargs(cfg),
        cfg.nstep,
        _ -> control_cfg.x0,
    )

    return OP.CenteredAbstractionGenerator{
        typeof(problem),
        typeof(gen_cfg),
        Any,
        Any,
    }(
        problem,
        gen_cfg,
        nothing,
        nothing,
        false,
        0.0,
    )
end

function build_certifier(problem, system_cfg, cfg::MarcheArriereConfig)
    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
    )

    backend = build_backend(; verbose = cfg.verbose)
    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        backend,
        opts,
    )

    params = system_cfg.params
    symbolic_builder = function (prob, candidate, certifier_cfg)
        o = certifier_cfg.options
        return AV.symbolic_system(
            prob.system.X;
            _U_ = prob.system.U,
            params = params,
            Ts = candidate.Ts,
            ΔX = o.ΔX,
            ΔU = o.ΔU,
            ΔW = o.ΔW,
            rk4_num_substeps = o.symbolic_rk4_substeps,
        )
    end

    return SC.EllipsoidalBackwardCertifier{
        typeof(problem),
        Any,
        typeof(cert_cfg),
        Any,
        typeof(symbolic_builder),
    }(
        nothing,
        nothing,
        cert_cfg,
        nothing,
        false,
        0.0,
        symbolic_builder,
    )
end

"""
Run generator + certification, forcing unwrapped certification trajectory when needed.
"""
function solve_with_unwrapped_certification!(solver::OP.CertifiedPipelineSolver, cfg::MarcheArriereConfig)
    @assert solver.problem !== nothing "Call set_problem!(solver, problem) first."

    candidate_for_cert = nothing
    certification_result = nothing
    success = false

    solve_time = @elapsed begin
        OP.generate!(solver.generator)
        candidate_raw = OP.get_trajectory(solver.generator)
        candidate_raw === nothing && error("Generator returned no candidate trajectory.")

        candidate_for_cert, was_unwrapped, wrap_jumps = candidate_for_periodic_certification(
            candidate_raw,
            cfg.periodic_dims,
            cfg.periodic_periods,
        )
         
        ## j'ai un probl!me avec la traj
        println("wad : ",was_unwrapped)
 
        if was_unwrapped
            println("Candidate unwrapped for certification: ", wrap_jumps, " periodic jumps fixed.")
        else
            println("Candidate already periodic-consistent for certification.")
        end

        SC.set_trajectory!(solver.certifier, candidate_for_cert)
        SC.certify!(solver.certifier)
        certification_result = SC.get_result(solver.certifier)
        certification_result === nothing && error("Certifier returned no result.")

        success = OP.get_success(solver.generator) && SC.get_success(solver.certifier)
    end

    solver.result = OP.CertifiedSolveResult(candidate_for_cert, certification_result, success, solve_time)
    return solver
end

function main(cfg::MarcheArriereConfig = MarcheArriereConfig())
    paths = output_paths(cfg)

    # ------------------------------------------------------------
    # 1) Build the different structure
    # ------------------------------------------------------------
    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)
    gen = build_generator(problem, system_cfg, control_cfg, cfg)
    cert = build_certifier(problem, system_cfg, cfg)

    solver = OP.CertifiedPipelineSolver{
        typeof(gen),
        typeof(cert),
        typeof(problem),
        Any,
    }(gen, cert, nothing, nothing)
    
    # ------------------------------------------------------------
    # 2) Build a candidate trajectory
    # ------------------------------------------------------------
    OP.set_problem!(solver, problem)
    solve_with_unwrapped_certification!(solver, cfg)
    result = OP.get_result(solver)

    result === nothing && error("Pipeline returned no result.")
    result.candidate === nothing && error("Pipeline produced no candidate trajectory.")

    save_state_space_plots!(
        paths.plots_dir,
        problem,
        result.candidate;
        cert_result = result.certification,
        show_ellipsoids = true,
        unwrap_angles = false,
        wrap_angles = true,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "marche_arriere (x,y)",
        title34 = "marche_arriere (theta,phi)",
    )

    gif_path = joinpath(paths.animations_dir, "rollout.gif")
    if cfg.plot_gif
        plot_articulated_vehicle!(
            AV,
            problem.system,
            system_cfg.params,
            result.candidate.x_traj,
            result.candidate.u_traj;
            giffile = gif_path,
            dt = result.candidate.Ts,
            title = "marche arriere - pipeline certifie",
        )
    end

    println("generator_success = ", OP.get_success(gen))
    println("certifier_success = ", SC.get_success(cert))
    println("pipeline_success = ", OP.get_success(solver))
    println("candidate_horizon = ", OP.horizon(result.candidate))
    if result.certification !== nothing
        println("cert_steps = ", length(result.certification.steps))
        println("failed_k = ", result.certification.failed_k)
    end
    println("output_root = ", paths.root)
    println("plots_dir = ", paths.plots_dir)
    cfg.plot_gif && println("gif = ", gif_path)

    return (
        ;
        solver,
        result,
        problem,
        config = cfg,
        outputs = paths,
        gif = cfg.plot_gif ? gif_path : nothing,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
