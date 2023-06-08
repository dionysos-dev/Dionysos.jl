export SCOTS

module SCOTS

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

using JuMP

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    problem::Union{Nothing, PR.ProblemType}
    symmodel::Union{Nothing, SY.SymbolicModelList}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2,NTuple{2,Int}}}
    controller
    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing
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

function get_abstract_system(optimizer)
    return optimizer.symmodel
end

function get_abstract_problem(optimizer)
    return optimizer.abstract_problem
end

function get_concrete_controller(optimizer)
    return optimizer.controller
end

function build_abstraction(
    system,
    state_grid::DO.Grid,
    input_grid::DO.Grid,
)
    Xfull = DO.DomainList(state_grid)
    DO.add_set!(Xfull, system.X, DO.INNER)
    Ufull = DO.DomainList(input_grid)
    DO.add_set!(Ufull, system.U, DO.CENTER)
    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    @time SY.compute_symmodel_from_controlsystem!(
        symmodel,
        system.f,
    )
    return symmodel
end


function build_abstract_problem(
    problem::PR.OptimalControlProblem,
    symmodel::SY.SymbolicModelList,
)
    state_grid = symmodel.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xtarget, symmodel.Xdom, problem.target_set, DO.INNER)
    return PR.OptimalControlProblem(
        symmodel,
        Xinit,
        Xtarget,
        problem.state_cost, # TODO this is the continuous cost, not the abstraction
        problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    problem::PR.SafetyProblem,
    symmodel::SY.SymbolicModelList,
)
    state_grid = symmodel.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, DO.OUTER)
    Xsafe = DO.DomainList(state_grid)
    DO.add_subset!(Xsafe, symmodel.Xdom, problem.safe_set, DO.INNER)
    return PR.SafetyProblem(
        symmodel,
        Xinit,
        Xsafe,
        problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function _compute_controller!(controller, discrete_problem::PR.OptimalControlProblem)
    init_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.initial_set)]
    target_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.target_set)]
    CO.compute_controller_reach!(
        controller,
        discrete_problem.system.autom,
        init_list,
        target_list,
    )
    return
end

function _compute_controller!(controller, discrete_problem::PR.SafetyProblem)
    init_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.initial_set)]
    safe_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.safe_set)]
    CO.compute_controller_safe!(
        controller,
        discrete_problem.system.autom,
        init_list,
        safe_list,
    )
    return
end

function MOI.optimize!(optimizer::Optimizer)
    # design the abstraction
    symmodel = build_abstraction(optimizer.problem.system, optimizer.state_grid, optimizer.input_grid)
    optimizer.symmodel = symmodel
    # design the abstract problem
    abstract_problem = build_abstract_problem(optimizer.problem, symmodel)
    optimizer.abstract_problem = abstract_problem
    # solve the abstract problem
    optimizer.abstract_controller = CO.NewControllerList()
    @time _compute_controller!(
        optimizer.abstract_controller,
        abstract_problem,
    )
    # refine the abstract controllerto the concrete controller
    function controller(x;param=false)
        xpos = DO.get_pos_by_coord(symmodel.Xdom.grid, x)
        if !(xpos ∈ symmodel.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(symmodel, xpos)
        symbollist = UT.fix_and_eliminate_first(optimizer.abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        upos = SY.get_upos_by_symbol(symmodel, symbol)
        u = DO.get_coord_by_pos(symmodel.Udom.grid, upos)
        return u
    end
    optimizer.controller = controller
    return 
end

end