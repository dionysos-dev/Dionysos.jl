export Abstraction

module Abstraction

using JuMP

import Dionysos

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, Dionysos.Domain.Grid}
    input_grid::Union{Nothing, Dionysos.Domain.Grid}
    problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem}
    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
        )
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end

function build_abstraction(
    system,
    state_grid::Dionysos.Domain.Grid,
    input_grid::Dionysos.Domain.Grid,
)
    Xfull = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_set!(Xfull, system.X, Dionysos.Domain.OUTER)
    Ufull = Dionysos.Domain.DomainList(input_grid)
    Dionysos.Domain.add_set!(Ufull, system.U, Dionysos.Domain.OUTER)
    symmodel = Dionysos.Symbolic.NewSymbolicModelListList(Xfull, Ufull)
    @time Dionysos.Symbolic.compute_symmodel_from_controlsystem!(
        symmodel,
        system,
    )
    return symmodel
end

function build_abstraction(
    problem::Dionysos.Problem.OptimalControlProblem,
    state_grid::Dionysos.Domain.Grid,
    input_grid::Dionysos.Domain.Grid,
)
    symmodel = build_abstraction(problem.system, state_grid, input_grid)
    Xinit = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, Dionysos.Problem.OUTER)
    Xtarget = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xinit, symmodel.Xdom, problem.target_set, Dionysos.Problem.OUTER)
    return Dionysos.Problem.OptimalControlProblem(
        symmodel,
        Xinit,
        Xtarget,
        problem.state_cost, # TODionysos.Problem this is the continuous cost, not the abstraction
        problem.transition_cost, # TODionysos.Problem this is the continuous cost, not the abstraction
        problem.time, # TODionysos.Problem this is the continuous time, not the number of transition
    )
end

function MOI.optimize!(optimizer::Optimizer{T}) where {T}
    discrete_problem = build_abstraction(optimizer.problem, optimizer.state_grid, optimizer.input_grid)
    optimizer.controller = Dionysos.Control.NewControllerList()
    @time Dionysos.Control.compute_controller_reach!(
        optimizer.controller,
        discrete_problem.system.autom,
        discrete_problem.initial_set,
        discrete_problem.target_set,
    )
    return
end

end