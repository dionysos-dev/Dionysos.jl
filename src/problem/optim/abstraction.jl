export Abstraction

module Abstraction

using JuMP

using ..Problem
using ...Dionysos.Control
using ...Dionysos.Domain

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, Domain.Grid}
    input_grid::Union{Nothing, Domain.Grid}
    problem::Union{Nothing, Problem.OptimalControlProblem}
    function Optimizer{T}() where {T}
        return new{T}(
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

function build_abstraction(system, state_grid::Domain.Grid, input_grid::Domain.Grid)
    Xfull = DO.DomainList(state_grid)
    DO.add_set!(Xfull, system.X, DO.OUTER)
    Ugrid = DO.Domain.GridFree(u0, h)
    Ufull = DO.DomainList(input_grid)
    DO.add_set!(Ufull, _U_, DO.OUTER)
    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    @time SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
    return symmodel
end

function build_abstraction(problem::Problem.OptimalControlProblem, state_grid::Domain.Grid, input_grid::Domain.Grid)
    symmodel = build_abstraction(problem.system, state_grid, input_grid)
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.target_set, DO.OUTER)
    return Problem.OptimalControlProblem(
        symmodel,
        Xinit,
        Xtarget,
        problem.state_cost, # TODO this is the continuous cost, not the abstraction
        problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function MOI.optimize!(optimizer::Optimizer{T}) where {T}
    discrete_problem = build_abstraction(optimizer.problem, optimizer.state_grid, optimizer.input_grid)
    optimizer.controller = CO.NewControllerList()
    @time CO.compute_controller_reach!(
        optimizer.controller,
        discrete_problem.system.autom,
        discrete_problem.initial_set,
        discrete_problem.target_set,
    )
    return
end

end