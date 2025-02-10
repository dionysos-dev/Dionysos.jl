export UniformGridAbstraction

module UniformGridAbstraction

import StaticArrays

import MathematicalSystems

import Dionysos

using JuMP

include("empty_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")

"""
    Optimizer{T} <: MOI.AbstractOptimizer

    A solver implementing the classical abstraction method (e.g., used in SCOTS), where the entire domain is partitioned into hyper-rectangular cells, independent of the specific control task. 

    The optimizer is structured into modular sub-solvers, each dedicated to a specific problem type. It ensures that abstraction is computed before solving the control problem, maintaining modularity and extendability.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    abstraction_solver::Union{Nothing, OptimizerEmptyProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Any
    solve_time_sec::T

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, 0.0)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.abstraction_solver === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    param_symbol = Symbol(param.name)

    if param_symbol == :concrete_problem
        # Assign appropriate control solver
        if isa(value, Dionysos.Problem.EmptyProblem)
            model.abstraction_solver = OptimizerEmptyProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                value,
            )
            model.control_solver = nothing  # Pas de solveur de contrôle
        elseif isa(value, Dionysos.Problem.OptimalControlProblem)
            model.control_solver = OptimizerOptimalControlProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        elseif isa(value, Dionysos.Problem.SafetyProblem)
            model.control_solver = OptimizerSafetyProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        else
            error("Unsupported problem type: $(typeof(value))")
        end
        # Instantiate an abstraction_solver if it has not already been created
        if model.abstraction_solver === nothing
            model.abstraction_solver = OptimizerEmptyProblem()
            empty_problem = Dionysos.Problem.EmptyProblem(value.system, value.system.X)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                empty_problem,
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

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

NewControllerList() = Dionysos.Utils.SortedTupleSet{2, NTuple{2, Int}}()

"""
    solve_concrete_problem(abstract_system, abstract_controller)

Solves the concrete problem based on the `abstract_system` and `abstract_controller`.
"""
function solve_concrete_problem(abstract_system, abstract_controller)
    function concrete_controller(x; param = false)
        # Getting the position of the state in the abstract system
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain: $x")
            return nothing
        end

        # Getting the corresponding abstract state
        source = Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
        symbollist = Dionysos.Utils.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state: $x")
            return nothing
        end

        # Choosing a random symbol or the first one
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end

        # Getting and return the control points
        upos = Dionysos.Symbolic.get_upos_by_symbol(abstract_system, symbol)
        u = Dionysos.Domain.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u
    end
end

function is_abstraction_computed(optimizer::Optimizer)
    return optimizer.abstraction_solver !== nothing &&
           optimizer.abstraction_solver.abstract_system !== nothing
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Ensure the concrete problem is defined
    if optimizer.abstraction_solver === nothing
        error("The concrete problem is not defined.")
    end

    # Compute abstraction if not already done
    if !is_abstraction_computed(optimizer)
        MOI.optimize!(optimizer.abstraction_solver)
    end

    # If there's a control solver, optimize it
    if optimizer.control_solver !== nothing
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
        optimizer.concrete_controller = solve_concrete_problem(
            optimizer.abstraction_solver.abstract_system,
            abstract_controller,
        )
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

end
