module UniformEllipsoidAbstraction

import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems

import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA
import Polyhedra
import LazySets
import MathematicalSystems as MS

using ProgressMeter
using JuMP
import MathOptInterface as MOI

export Optimizer

include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")

"""
    Optimizer{T} <: Dionysos.Optim.AbstractionControlOptimizer

Abstraction-based solver using a uniform partition of ellipsoidal cells: it builds the ellipsoidal
abstraction, synthesizes an abstract controller (forwarding the transition-cost closure), and
concretizes it into a `RefinedStaticController`. Follows the shared
[`AbstractionControlOptimizer`](@ref Dionysos.Optim.AbstractionControlOptimizer) pipeline.
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

function OP.configure_control_solver!(model::Optimizer, control_solver, abstract_system)
    transitionCost = model.abstraction_solver.transitionCost
    costfun = (q::Int, s::Int) -> get(transitionCost, (q, s), Inf)
    MOI.set(control_solver, MOI.RawOptimizerAttribute("abstract_transition_cost"), costfun)
    return
end

function OP.build_concrete_controller(
    model::Optimizer,
    abstract_system,
    abstract_controller,
)
    abstract_value_function =
        MOI.get(model.control_solver, MOI.RawOptimizerAttribute("abstract_value_function"))
    transitionCont =
        MOI.get(model.abstraction_solver, MOI.RawOptimizerAttribute("transitionCont"))
    return RefinedStaticController(
        model.abstraction_solver.abstract_system,
        abstract_controller,
        transitionCont,
        abstract_value_function,
    )
end

struct RefinedStaticController{AS, AC, TC, VF} <: ST.AbstractContinuousController
    abstract_system::AS
    abstract_controller::AC
    transition_controllers::TC
    abstract_value_function::VF
    randomize::Bool
end

function RefinedStaticController(
    abstract_system,
    abstract_controller::ST.AbstractDiscreteController,
    transitionCont::Dict,
    abstract_value_function::Function,
)
    return RefinedStaticController(
        abstract_system,
        abstract_controller,
        transitionCont,
        abstract_value_function,
        false,
    )
end

ST.controller_kind(::RefinedStaticController) = ST.StaticKind()
ST.domain(ctrl::RefinedStaticController) = ctrl.abstract_system

function ST.is_defined(ctrl::RefinedStaticController, q, x)
    return pick_best_refined_transition(ctrl, x) !== nothing
end

function ST.output_control(ctrl::RefinedStaticController, q, x)
    tr = pick_best_refined_transition(ctrl, x)
    tr === nothing && return nothing
    return MS.apply(tr.cont, tr.x)
end

function pick_best_refined_transition(ctrl::RefinedStaticController, x)
    qs = SY.get_abstract_states(ctrl.abstract_system, x)
    isempty(qs) && return nothing

    best_q = nothing
    best_to = nothing
    best_cont = nothing
    best_val = Inf

    for q in qs
        ST.is_defined(ctrl.abstract_controller, nothing, q) || continue

        to_list = ctrl.abstract_controller.controller_map(q)
        to_list === nothing && continue
        isempty(to_list) && continue

        to = ctrl.randomize ? rand(to_list) : first(to_list)

        key = (q, to)
        haskey(ctrl.transition_controllers, key) || continue

        v = ctrl.abstract_value_function(q)
        if v < best_val
            best_val = v
            best_q = q
            best_to = to
            best_cont = ctrl.transition_controllers[key]
        end
    end

    best_q === nothing && return nothing

    return (x = x, from = best_q, to = best_to, cont = best_cont, val = best_val)
end

end # module
