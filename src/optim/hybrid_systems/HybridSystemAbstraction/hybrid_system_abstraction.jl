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

include("hybrid_symbolic_builder.jl")
include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")

"""
    Optimizer{T} <: Dionysos.Optim.AbstractionControlOptimizer

Abstraction-based solver for timed hybrid systems: it builds the
[`HybridSymbolicModel`](@ref Dionysos.Symbolic.HybridSymbolicModel) abstraction, solves the
abstract control problem, and concretizes the controller into a `HybridQuantizedStaticController`.
Follows the shared [`AbstractionControlOptimizer`](@ref Dionysos.Optim.AbstractionControlOptimizer) pipeline.
"""
mutable struct Optimizer{T} <: OP.AbstractionControlOptimizer
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

OP.default_abstraction_solver(::Optimizer) = OptimizerAlternatingSimulationProblem()

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

function OP.build_concrete_controller(::Optimizer, abstract_system, abstract_controller)
    return HybridQuantizedStaticController(abstract_system, abstract_controller)
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
    abs_q <= 0 && return false
    return ST.is_defined(ctrl.abstract_controller, nothing, abs_q)
end

function ST.output_control(ctrl::HybridQuantizedStaticController, q, aug_state)
    abs_q = SY.get_abstract_state(ctrl.abstract_system, aug_state)
    abs_q <= 0 && return nothing

    u_abs = ST.output_control(ctrl.abstract_controller, nothing, abs_q)
    u_abs === nothing && return nothing

    k = aug_state[end]
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

# Locate the hybrid transition encoded by a "SWITCH src -> tgt" label.
function _find_switch_transition(hs::HybridSystem, u)
    m = match(r"SWITCH (\d+) -> (\d+)", String(u))
    m === nothing && error("Unrecognized SWITCH label format: $u")
    source_mode = parse(Int, m.captures[1])
    target_mode = parse(Int, m.captures[2])

    transitions = collect(HybridSystems.transitions(hs.automaton))
    transition_id = findfirst(
        tr ->
            HybridSystems.source(hs.automaton, tr) == source_mode &&
            HybridSystems.target(hs.automaton, tr) == target_mode,
        transitions,
    )
    transition_id === nothing && error("Transition not found for $u")
    return transitions[transition_id], target_mode
end

# Advance one closed-loop step. The augmented state is `(x, t, mode)` for a
# clock-lifted mode or `(x, mode)` for a time-free mode; the time slot is advanced
# (or reset by the guard's reset map) only in the timed case.
function get_next_aug_state(hs::HybridSystem, aug_state, u, time_is_active, tstep, map_sys)
    x = aug_state[1]
    k = aug_state[end]
    timed = length(aug_state) == 3

    if isa(u, AbstractString) && occursin("SWITCH", u)
        transition, target_mode = _find_switch_transition(hs, u)
        reset_map = HybridSystems.resetmap(hs, transition)
        if timed
            t = aug_state[2]
            reset_result = MathematicalSystems.apply(reset_map, vcat(x, t))
            next_x = reset_result[1:(end - 1)]
            next_t = round(reset_result[end]; digits = 10)
            return (next_x, next_t, target_mode)
        else
            return (MathematicalSystems.apply(reset_map, x), target_mode)
        end
    else
        next_x = map_sys(x, u, tstep)
        if timed
            t = aug_state[2]
            next_t = round(time_is_active ? t + tstep : 0.0; digits = 10)
            return (next_x, next_t, k)
        else
            return (next_x, k)
        end
    end
end

# Physical dynamics `f` of a mode: `systems[1]` for a `VectorContinuousSystem`
# (paired with a clock), or the mode system itself when it is time-free.
_mode_dynamics(mode_system) =
    mode_system isa ST.VectorContinuousSystem ? mode_system.systems[1].f : mode_system.f

# Whether a mode's clock evolves (`ẋ = 1`); a time-free mode has no clock.
_mode_time_active(mode_system) =
    mode_system isa ST.VectorContinuousSystem && ([1.0;;] == mode_system.systems[2].A)

# Build the augmented-state step map `(aug, u) -> next_aug` for a hybrid system,
# hoisting the per-mode dynamics / clock-activity / time-step lookups once.
function _hybrid_step_map(hs::HybridSystem, tsteps)
    nmodes = HybridSystems.nmodes(hs.automaton)
    mode_systems = [HybridSystems.mode(hs, k) for k in 1:nmodes]
    maps_sys = [ST.simulate_control_map(_mode_dynamics(ms)) for ms in mode_systems]
    times_is_active = [_mode_time_active(ms) for ms in mode_systems]

    return function (aug_state, u)
        k = aug_state[end]
        return get_next_aug_state(
            hs,
            aug_state,
            u,
            times_is_active[k],
            tsteps[k],
            maps_sys[k],
        )
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
    # Reuse the generic closed-loop engine: the "state" is the augmented
    # `(x[, t], mode)`, stepped by `_hybrid_step_map`, and inputs are
    # heterogeneous (continuous inputs mixed with switching labels), hence
    # `input_type = Any`.
    traj = ST.get_closed_loop_trajectory(
        hs,
        controller,
        aug_state_0,
        nstep;
        stopping = stopping,
        f_map_override = _hybrid_step_map(hs, tsteps),
        input_type = Any,
    )
    return ST.states(traj), ST.inputs(traj)
end

end # module
