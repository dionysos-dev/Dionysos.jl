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
const SC = OP.SymbolicCertifier

# Compat refactor: certains problemes appellent encore SY.format_* alors que
# les helpers sont maintenant exposes dans Utils.
if !isdefined(DI.Symbolic, :format_input_set)
    @eval DI.Symbolic format_input_set(args...) = $(UT).format_input_set(args...)
end
if !isdefined(DI.Symbolic, :format_noise_set)
    @eval DI.Symbolic format_noise_set(args...) = $(UT).format_noise_set(args...)
end

# Charge le probleme articule exactement depuis Dionysos.
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

# Helpers de plot utilises dans Permis_E.
include(joinpath(@__DIR__, "helpers", "helpers_plots.jl"))

# (c'est pas très joli pour le moment c'est pratique, je voudrais améliorer l'ergo, comment faire  (?))
""" 
Configuration compacte du pipeline:
1) generation d'une trajectoire candidate,
2) certification ellipsoïdale backward,
3) export des plots (pdf + gif).
"""
Base.@kwdef struct SolverConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 700

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
        IA.interval(-0.2, 0.2),
        IA.interval(-0.2, 0.2),
    )
    ΔW::IA.IntervalBox{1, Float64} = IA.IntervalBox(IA.interval(0.0, 0.0), 1)

    outdir::String = joinpath(@__DIR__, "outputs", "articulated_pipeline")
    plot_gif::Bool = true
    verbose::Bool = true
end

periodicity_kwargs(cfg::SolverConfig) = (;
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

function build_backend(; verbose::Bool = false) # utilise mosek ou clarabel
    return optimizer_with_attributes(
        MosekTools.Optimizer,
        MOI.Silent() => !verbose,
        # 10 = detailed log, 0 = silent.
        MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => (verbose ? 10 : 0),
        MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
    )
end

function build_concrete_system()
    x_domain = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi),
        SVector(10.0, 9.0, pi, pi),
    )
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(50.0))

    obstacles_xy = [
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 6.0), SVector(10.0, 9.0)),
    ]
    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = pi / 4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max))

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)
    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_mapping()
    inputs_delta = [
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
        SVector(9.0, 5.0, pi - 5 * (pi / 180), -5 * (pi / 180)),
        SVector(10.0, 6.0, pi + 5 * (pi / 180), 5 * (pi / 180)),
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

function build_generator(problem, system_cfg, control_cfg, cfg::SolverConfig)
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

function build_certifier(problem, system_cfg, cfg::SolverConfig)
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

function _count_periodic_wrap_jumps(state_list, periodic_dims, periodic_periods)
    length(state_list) <= 1 && return 0
    njump = 0
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        threshold = 0.5 * p
        for k in 2:length(state_list)
            Δ = Float64(state_list[k][d] - state_list[k - 1][d])
            abs(Δ) > threshold && (njump += 1)
        end
    end
    return njump
end

function _augment_metadata_for_unwrap(metadata, was_unwrapped::Bool, wrap_jumps::Int)
    if metadata isa NamedTuple
        return merge(
            metadata,
            (;
                lmi_unwrapped = was_unwrapped,
                lmi_wrap_jumps = wrap_jumps,
            ),
        )
    end
    return (;
        original_metadata = metadata,
        lmi_unwrapped = was_unwrapped,
        lmi_wrap_jumps = wrap_jumps,
    )
end

function candidate_for_certification(candidate::OP.CandidateTrajectory, cfg::SolverConfig)
    xs = collect(ST.enum_elems(candidate.x_traj))
    wrap_jumps = _count_periodic_wrap_jumps(xs, cfg.periodic_dims, cfg.periodic_periods)

    if wrap_jumps == 0
        return candidate, false, wrap_jumps
    end

    xs_unwrapped = unwrap_periodic_state_list(xs, cfg.periodic_dims, cfg.periodic_periods)
    cand_unwrapped = OP.CandidateTrajectory(
        ST.Trajectory(xs_unwrapped),
        candidate.u_traj;
        Ts = candidate.Ts,
        source = candidate.source,
        metadata = _augment_metadata_for_unwrap(candidate.metadata, true, wrap_jumps),
    )
    return cand_unwrapped, true, wrap_jumps
end

function solve_with_unwrapped_certification!(
    solver::OP.CertifiedPipelineSolver,
    cfg::SolverConfig,
)
    @assert solver.problem !== nothing "Call set_problem!(solver, problem) first."

    cand_for_cert = nothing
    certres = nothing
    ok = false

    t = @elapsed begin
        OP.generate!(solver.generator)
        cand_raw = OP.get_trajectory(solver.generator)
        cand_raw === nothing && error("Aucune candidate generee.")

        cand_for_cert, was_unwrapped, wrap_jumps = candidate_for_certification(cand_raw, cfg)
        if was_unwrapped
            println("Candidate unwrap pour certif: ", wrap_jumps, " saut(s) periodiques corriges.")
        else
            println("Candidate deja unwrap-coherente pour certif (aucun saut periodique detecte).")
        end

        SC.set_trajectory!(solver.certifier, cand_for_cert)
        SC.certify!(solver.certifier)
        certres = SC.get_result(solver.certifier)
        certres === nothing && error("Le certifier n'a retourne aucun resultat.")

        ok = OP.get_success(solver.generator) && SC.get_success(solver.certifier)
    end

    solver.result = OP.CertifiedSolveResult(cand_for_cert, certres, ok, t)
    return solver
end

"""
Pipeline principal:
- construit le probleme articulated vehicle,
- genere une candidate (CenteredAbstractionGenerator),
- certifie (EllipsoidalBackwardCertifier via CertifiedPipelineSolver),
- sauve les plots.
"""
function main(cfg::SolverConfig = SolverConfig())
    mkpath(cfg.outdir)

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

    OP.set_problem!(solver, problem)
    solve_with_unwrapped_certification!(solver, cfg)
    res = OP.get_result(solver)
    res === nothing && error("Le pipeline n'a retourne aucun resultat.")
    res.candidate === nothing && error("Aucune trajectoire candidate n'a ete produite.")

    save_state_space_plots!(
        cfg.outdir,
        problem,
        res.candidate;
        cert_result = res.certification,
        show_ellipsoids = true,
        unwrap_angles = false,
        wrap_angles = true,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        title12 = "articulated_vehicle (x,y)",
        title34 = "articulated_vehicle (theta,phi)",
    )

    gifpath = joinpath(cfg.outdir, "rollout.gif")
    if cfg.plot_gif
        plot_articulated_vehicle!(
            AV,
            problem.system,
            system_cfg.params,
            res.candidate.x_traj,
            res.candidate.u_traj;
            giffile = gifpath,
            dt = res.candidate.Ts,
            title = "Articulated vehicle - pipeline certifie",
        )
    end

    println("generator_success = ", OP.get_success(gen))
    println("certifier_success = ", SC.get_success(cert))
    println("pipeline_success = ", OP.get_success(solver))
    println("candidate_horizon = ", OP.horizon(res.candidate))
    if res.certification !== nothing
        println("cert_steps = ", length(res.certification.steps))
        println("failed_k = ", res.certification.failed_k)
    end
    println("plots_dir = ", cfg.outdir)
    cfg.plot_gif && println("gif = ", gifpath)

    return (; solver, result = res, problem, config = cfg, outdir = cfg.outdir, gif = gifpath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
