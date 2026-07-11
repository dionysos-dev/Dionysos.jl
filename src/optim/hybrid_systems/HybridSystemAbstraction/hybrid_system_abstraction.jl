export HybridSystemAbstraction

module HybridSystemAbstraction
using MathOptInterface
using HybridSystems
import MathematicalSystems
MS = MathematicalSystems
using StaticArrays: SVector
using Dionysos
import LazySets

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems

export Optimizer

include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")

mutable struct Optimizer{T} <: OP.CompositeDionysosOptimizer
    abstraction_solver::Union{Nothing, OptimizerAlternatingSimulationProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, ST.AbstractContinuousController}
    solve_time_sec::T
    print_level::Int

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, 0.0, 1)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.abstraction_solver === nothing

OP.sub_solvers(model::Optimizer) = (model.abstraction_solver, model.control_solver)

function OP.ensure_sub_solvers!(model::Optimizer)
    if model.abstraction_solver === nothing
        model.abstraction_solver = OptimizerAlternatingSimulationProblem()
    end
    return
end

"""
    control_solver_for(problem)

Return a fresh control sub-solver matching the type of the concrete `problem`. Supporting a new
`ProblemType` in this solver family means adding one method here.
"""
control_solver_for(::PR.OptimalControlProblem) = OptimizerOptimalControlProblem()
control_solver_for(::PR.SafetyProblem) = OptimizerSafetyProblem()
control_solver_for(problem) = error("Unsupported problem type: $(typeof(problem))")

function OP.set_concrete_problem!(
    model::Optimizer,
    problem::PR.AlternatingSimulationProblem,
)
    model.abstraction_solver = OptimizerAlternatingSimulationProblem()
    model.abstraction_solver.alternating_simulation_problem = problem
    model.control_solver = nothing
    return
end

function OP.set_concrete_problem!(model::Optimizer, problem)
    model.control_solver = control_solver_for(problem)
    model.control_solver.concrete_problem = problem
    if model.abstraction_solver.alternating_simulation_problem === nothing
        model.abstraction_solver.alternating_simulation_problem =
            PR.AlternatingSimulationProblem(problem.system, nothing)
    end
    return
end

function is_abstraction_computed(optimizer::Optimizer)
    return optimizer.abstraction_solver !== nothing &&
           optimizer.abstraction_solver.abstract_system !== nothing
end

function reset!(optimizer::Optimizer)
    optimizer.concrete_controller = nothing
    optimizer.solve_time_sec = 0.0
    if optimizer.control_solver !== nothing
        reset!(optimizer.control_solver)
    end
    if optimizer.abstraction_solver !== nothing
        reset!(optimizer.abstraction_solver)
    end
    return optimizer
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Ensure the concrete problem is defined
    if optimizer.abstraction_solver === nothing
        error("The concrete problem is not defined.")
    end

    # Compute abstraction if not already done
    if !is_abstraction_computed(optimizer)
        MOI.set(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        MOI.optimize!(optimizer.abstraction_solver)
    end

    # If there's a control solver, optimize it
    if optimizer.control_solver !== nothing
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        abstract_system = MOI.get(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
        )
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            abstract_system,
        )
        MOI.optimize!(optimizer.control_solver)
        abstract_controller = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_controller"),
        )
        optimizer.concrete_controller =
            HybridQuantizedStaticController(abstract_system, abstract_controller)
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

# ================================================================
# Concrete controller synthesis
# ================================================================

struct HybridQuantizedStaticController{AS, AC} <: ST.AbstractContinuousController
    abstract_system::AS
    abstract_controller::AC
end
ST.controller_kind(::HybridQuantizedStaticController) = ST.StaticKind()
ST.domain(ctrl::HybridQuantizedStaticController) = ctrl.abstract_system

function ST.is_defined(ctrl::HybridQuantizedStaticController, q, aug_state)
    abs_q = SY.get_abstract_state(ctrl.abstract_system, aug_state)
    abs_q === nothing && return false
    return ST.is_defined(ctrl.abstract_controller, nothing, abs_q)
end

function ST.output_control(ctrl::HybridQuantizedStaticController, q, aug_state)
    abs_q = SY.get_abstract_state(ctrl.abstract_system, aug_state)
    abs_q === nothing && return nothing

    u_abs = ST.output_control(ctrl.abstract_controller, nothing, abs_q)
    u_abs === nothing && return nothing

    (_, _, k) = aug_state
    if SY.is_switching_input(ctrl.abstract_system.input_mapping, u_abs)
        transition_id = ctrl.abstract_system.input_mapping.global_to_switching[u_abs]
        return ctrl.abstract_system.input_mapping.switch_labels[transition_id]
    else
        return SY.get_concrete_input(ctrl.abstract_system, u_abs, k)
    end
end

# ================================================================
# Closed-loop simulation utilities
# ================================================================

function get_next_aug_state(hs::HybridSystem, aug_state, u, time_is_active, tstep, map_sys)
    (x, t, k) = aug_state

    if isa(u, AbstractString) && occursin("SWITCH", u)
        m = match(r"SWITCH (\d+) -> (\d+)", String(u))
        if m === nothing
            error("Unrecognized SWITCH label format: $u")
        end
        source_mode = parse(Int, m.captures[1])
        target_mode = parse(Int, m.captures[2])

        transitions = collect(HybridSystems.transitions(hs.automaton))
        transition_id = findfirst(
            tr ->
                HybridSystems.source(hs.automaton, tr) == source_mode &&
                HybridSystems.target(hs.automaton, tr) == target_mode,
            transitions,
        )
        if transition_id === nothing
            error("Transition not found for $u")
        end
        transition = transitions[transition_id]

        reset_map = HybridSystems.resetmap(hs, transition)
        augmented_source_state = vcat(x, t)
        reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
        next_x = reset_result[1:(end - 1)]
        next_t = reset_result[end]
        next_k = target_mode
        next_t = round(next_t; digits = 10)
        return (next_x, next_t, next_k)
    else
        next_t = time_is_active ? t + tstep : 0.0
        next_t = round(next_t; digits = 10)
        next_x = map_sys(x, u, tstep)
        return (next_x, next_t, k)
    end
end

function get_closed_loop_trajectory(
    hs::HybridSystem,
    controller::ST.AbstractController,
    tsteps,
    aug_state_0,
    nstep;
    stopping = x -> false,
)
    q = ST.initial_state(controller)
    aug_state_traj, u_traj = [aug_state_0], Any[]
    aug_state = aug_state_0

    nmodes = HybridSystems.nmodes(hs.automaton)
    dynamics = [HybridSystems.mode(hs, k).systems[1].f for k in 1:nmodes]
    maps_sys = [ST.simulate_control_map(dynamics[i]) for i in 1:nmodes]
    times_is_active =
        [([1.0;;] == HybridSystems.mode(hs, k).systems[2].A) for k in 1:nmodes]

    for _ in 1:nstep
        stopping(aug_state) && break

        u = ST.output_control(controller, q, aug_state)
        u === nothing && break

        (_, _, k) = aug_state
        next_aug_state =
            get_next_aug_state(hs, aug_state, u, times_is_active[k], tsteps[k], maps_sys[k])

        q = ST.update_state(controller, q, aug_state)
        aug_state = next_aug_state

        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end

    return aug_state_traj, u_traj
end

end # module
