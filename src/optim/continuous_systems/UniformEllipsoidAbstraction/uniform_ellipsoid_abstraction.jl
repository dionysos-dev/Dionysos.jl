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
import MathematicalSystems as MS

using ProgressMeter
using JuMP
import MathOptInterface as MOI

export Optimizer

include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
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

MOI.is_empty(opt::Optimizer) = opt.abstraction_solver === nothing

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    return model.print_level = value ? 0 : 1
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    param_symbol = Symbol(param.name)

    if model.abstraction_solver === nothing
        model.abstraction_solver = OptimizerAlternatingSimulationProblem()
    end

    if param_symbol == :abstract_system
        MOI.set(
            model.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            value,
        )
        return
    end

    if param_symbol == :concrete_problem
        # Assign appropriate control solver
        if isa(value, Dionysos.Problem.AlternatingSimulationProblem)
            model.abstraction_solver = OptimizerAlternatingSimulationProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                value,
            )
            model.control_solver = nothing  # No control solver
        elseif isa(value, Dionysos.Problem.OptimalControlProblem)
            model.control_solver = OptimizerOptimalControlProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        else
            error("Unsupported problem type: $(typeof(value))")
        end

        # Instantiate an abstraction_solver if it has not already been created
        if model.abstraction_solver.alternating_simulation_problem === nothing
            alternating_simulation_problem =
                Dionysos.Problem.AlternatingSimulationProblem(value.system, nothing)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                alternating_simulation_problem,
            )
        end
        return
    end

    # If the attribute exists in the main optimizer, set it there
    if hasfield(typeof(model), param_symbol)
        return setproperty!(model, param_symbol, value)
    end

    # Try setting it in the sub-solvers
    for solver in (model.abstraction_solver, model.control_solver)
        if solver !== nothing && hasfield(typeof(solver), param_symbol)
            return setproperty!(solver, param_symbol, value)
        end
    end

    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    param_symbol = Symbol(param.name)

    # First, check if the attribute exists directly in the main optimizer
    if hasfield(typeof(model), param_symbol)
        return getproperty(model, param_symbol)
    end

    # If not found, try to get it from the abstraction solver
    if model.abstraction_solver !== nothing &&
       hasfield(typeof(model.abstraction_solver), param_symbol)
        return getproperty(model.abstraction_solver, param_symbol)
    end

    # If not found, try to get it from the control solver
    if model.control_solver !== nothing &&
       hasfield(typeof(model.control_solver), param_symbol)
        return getproperty(model.control_solver, param_symbol)
    end

    # If the attribute is not recognized, raise an error
    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

MOI.get(opt::Optimizer, ::MOI.SolveTimeSec) = opt.solve_time_sec

function is_abstraction_computed(opt::Optimizer)
    return opt.abstraction_solver !== nothing &&
           opt.abstraction_solver.abstract_system !== nothing
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
        transitionCost = optimizer.abstraction_solver.transitionCost
        costfun = (q::Int, s::Int) -> get(transitionCost, (q, s), Inf)
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_transition_cost"),
            costfun,
        )

        MOI.optimize!(optimizer.control_solver)
        abstract_controller = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_controller"),
        )
        abstract_value_function = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_value_function"),
        )
        transitionCont = MOI.get(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("transitionCont"),
        )
        optimizer.concrete_controller = RefinedStaticController(
            optimizer.abstraction_solver.abstract_system,
            abstract_controller,
            transitionCont,
            abstract_value_function,
        )
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

struct RefinedStaticController{AS, AC, TC, VF, H} <: ST.AbstractContinuousController
    abstract_system::AS
    abstract_controller::AC
    transition_controllers::TC
    abstract_value_function::VF
    handle_out_of_domain::H
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
        (x, sys) -> nothing,
        false,
    )
end

ST.domain(ctrl::RefinedStaticController) = ctrl.abstract_system
ST.initial_state(::RefinedStaticController) = nothing
ST.update_state(::RefinedStaticController, x, y) = nothing

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

    if isempty(qs)
        xnew = ctrl.handle_out_of_domain(x, ctrl.abstract_system)
        xnew === nothing && return nothing
        qs = SY.get_abstract_states(ctrl.abstract_system, xnew)
        isempty(qs) && return nothing
        x = xnew
    end

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
