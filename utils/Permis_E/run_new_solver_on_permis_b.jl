using LinearAlgebra
import StaticArrays: SVector
import Dionysos
import MathOptInterface as MOI

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const HG = OP
const SC = OP.SymbolicCertifier
const CS = OP

# -------------------------------------------------------------------
# Compat layer (legacy Permis_B API -> refactored Dionysos API)
# -------------------------------------------------------------------

# Legacy namespace expected by utils/Permis_B/permis_B.jl.
if !isdefined(DI, :Domain)
    @eval DI module Domain
        const _MP = $(DI).Mapping
        const INNER = _MP.INNER
        const OUTER = _MP.OUTER
        const CENTER = _MP.CENTER

        GridFree(args...; kwargs...) = _MP.GridFree(args...; kwargs...)
        CustomList(inputs) = _MP.ListMapping(inputs)
        DomainList(grid) = _MP.ExplicitGridMapping(grid)
        add_set!(dom, set, incl_mode) = (_MP.add_set!(dom, set, incl_mode); dom)
    end
end

# Legacy names in articulated_vehicle symbolic builder.
if !isdefined(DI.Symbolic, :format_input_set)
    @eval DI.Symbolic format_input_set(args...) = $(UT).format_input_set(args...)
end
if !isdefined(DI.Symbolic, :format_noise_set)
    @eval DI.Symbolic format_noise_set(args...) = $(UT).format_noise_set(args...)
end

# Legacy MOI attributes used in Permis_B helper scripts.
const UGA = OP.Abstraction.UniformGridAbstraction

function MOI.set(
    model::UGA.Optimizer,
    attr::MOI.RawOptimizerAttribute,
    value::MP.AbstractMapping,
)
    if attr.name == "Udom"
        return MOI.set(model, MOI.RawOptimizerAttribute("UMapping"), value)
    end
    return invoke(
        MOI.set,
        Tuple{UGA.Optimizer, MOI.RawOptimizerAttribute, Any},
        model,
        attr,
        value,
    )
end

function MOI.set(
    model::UGA.Optimizer,
    attr::MOI.RawOptimizerAttribute,
    value::Bool,
)
    if attr.name == "use_periodic_domain"
        return MOI.set(model, MOI.RawOptimizerAttribute("use_periodic_mapping"), value)
    end
    return invoke(
        MOI.set,
        Tuple{UGA.Optimizer, MOI.RawOptimizerAttribute, Any},
        model,
        attr,
        value,
    )
end

module PBRef
include(joinpath(@__DIR__, "..", "Permis_B", "permis_B.jl"))
end

Base.@kwdef mutable struct PermisBCandidateGenerator
    sim_cfg::PBRef.SimulationConfig = PBRef.SimulationConfig()
    problem::Any = nothing
    candidate::Any = nothing
    success::Bool = false
    solve_time_sec::Float64 = 0.0
end

function OP.set_problem!(gen::PermisBCandidateGenerator, prob)
    gen.problem = prob
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    return gen
end

function OP.generate!(gen::PermisBCandidateGenerator)
    t = @elapsed begin
        run_result = PBRef.run_nominal_simulation(
            gen.sim_cfg;
            save_plots = false,
            show_animation = false,
            use_cache = true,
            force_recompute = false,
        )

        run_result_lmi = PBRef.prepare_run_result_for_lmi(
            run_result;
            enable_unwrap = true,
            periodic_dims = gen.sim_cfg.periodic_dims,
            periodic_periods = gen.sim_cfg.periodic_periods,
        )

        cand = OP.CandidateTrajectory(
            ST.Trajectory(run_result_lmi.state_list),
            ST.Trajectory(run_result_lmi.input_list);
            Ts = run_result.Δt,
            source = :permis_b_nominal,
            metadata = (; from_cache = hasproperty(run_result, :from_cache) ? run_result.from_cache : false),
        )

        gen.candidate = cand
        if gen.problem !== nothing && hasproperty(gen.problem, :target_set)
            wrap = ST.get_periodic_wrapper(
                gen.sim_cfg.periodic_dims,
                gen.sim_cfg.periodic_periods;
                start = gen.sim_cfg.periodic_start,
            )
            gen.success = wrap(last(run_result_lmi.state_list)) ∈ gen.problem.target_set
        else
            gen.success = true
        end
    end
    gen.solve_time_sec = t
    return gen
end

OP.get_trajectory(gen::PermisBCandidateGenerator) = gen.candidate
OP.get_success(gen::PermisBCandidateGenerator) = gen.success
OP.get_solve_time(gen::PermisBCandidateGenerator) = gen.solve_time_sec

# Bridge OP-level orchestration calls to SymbolicCertifier methods.
@eval OP begin
    set_problem!(
        cert::SymbolicCertifier.AbstractSymbolicCertifier,
        prob::Dionysos.Problem.ProblemType,
    ) = SymbolicCertifier.set_problem!(cert, prob)
    set_trajectory!(cert::SymbolicCertifier.AbstractSymbolicCertifier, cand) =
        SymbolicCertifier.set_trajectory!(cert, cand)
    certify!(cert::SymbolicCertifier.AbstractSymbolicCertifier) =
        SymbolicCertifier.certify!(cert)
    get_result(cert::SymbolicCertifier.AbstractSymbolicCertifier) =
        SymbolicCertifier.get_result(cert)
    get_success(cert::SymbolicCertifier.AbstractSymbolicCertifier) =
        SymbolicCertifier.get_success(cert)
    get_solve_time(cert::SymbolicCertifier.AbstractSymbolicCertifier) =
        SymbolicCertifier.get_solve_time(cert)
end

function main()
    sim_cfg = PBRef.SimulationConfig()
    lmi_cfg = PBRef.LMIConfig()

    system_cfg = PBRef.build_concrete_system()
    control_cfg = PBRef.build_control_problem()

    problem = PR.OptimalControlProblem(
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )

    gen = PermisBCandidateGenerator(; sim_cfg = sim_cfg)

    opts = (
        λ = lmi_cfg.λ,
        maxδx = lmi_cfg.maxδx,
        maxδu = lmi_cfg.maxδu,
        ΔX = lmi_cfg.ΔX,
        ΔU = lmi_cfg.ΔU,
        ΔW = lmi_cfg.ΔW,
        rayon_terminal = lmi_cfg.rayon_terminal,
        symbolic_rk4_substeps = lmi_cfg.symbolic_rk4_substeps,
    )

    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        lmi_cfg.sdp_opt,
        opts,
    )

    pb_params = system_cfg.params
    builder = function (prob, candidate, cfg)
        return PBRef.AV.symbolic_system(
            prob.system.X;
            _U_ = prob.system.U,
            params = pb_params,
            Ts = candidate.Ts,
            ΔX = cfg.options.ΔX,
            ΔU = cfg.options.ΔU,
            ΔW = cfg.options.ΔW,
            rk4_num_substeps = cfg.options.symbolic_rk4_substeps,
        )
    end

    cert = SC.EllipsoidalBackwardCertifier{
        typeof(problem),
        Any,
        typeof(cert_cfg),
        Any,
        typeof(builder),
    }(
        nothing,
        nothing,
        cert_cfg,
        nothing,
        false,
        0.0,
        builder,
    )

    solver = CS.CertifiedPipelineSolver{
        typeof(gen),
        typeof(cert),
        typeof(problem),
        Any,
    }(gen, cert, nothing, nothing)
    CS.set_problem!(solver, problem)
    CS.solve!(solver)
    res = CS.get_result(solver)

    println("generator success=", HG.get_success(gen), " time=", HG.get_solve_time(gen))
    println("certifier success=", SC.get_success(cert), " time=", SC.get_solve_time(cert))
    println("pipeline success=", CS.get_success(solver), " time=", CS.get_solve_time(solver))
    if res !== nothing
        println("n_steps_cert=", length(res.certification.steps))
        println("failed_k=", res.certification.failed_k)
    end

    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
