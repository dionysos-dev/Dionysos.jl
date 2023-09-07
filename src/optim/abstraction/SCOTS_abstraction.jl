export SCOTSAbstraction

module SCOTSAbstraction

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
    concrete_problem::Union{Nothing, PR.ProblemType}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Any
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    solve_time_sec::T
    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, nothing, nothing, nothing, nothing, 0.0)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstraction(concrete_system, state_grid::DO.Grid, input_grid::DO.Grid)
    Xfull = DO.DomainList(state_grid)
    DO.add_set!(Xfull, concrete_system.X, DO.INNER)
    Ufull = DO.DomainList(input_grid)
    DO.add_set!(Ufull, concrete_system.U, DO.CENTER)
    abstract_system = SY.NewSymbolicModelListList(Xfull, Ufull)
    @time SY.compute_symmodel_from_controlsystem!(abstract_system, concrete_system.f)
    return abstract_system
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xtarget, abstract_system.Xdom, concrete_problem.target_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    target_list =
        [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xtarget)]
    return PR.OptimalControlProblem(
        abstract_system,
        init_list,
        target_list,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xsafe = DO.DomainList(state_grid)
    DO.add_subset!(Xsafe, abstract_system.Xdom, concrete_problem.safe_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    safe_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xsafe)]
    return PR.SafetyProblem(
        abstract_system,
        init_list,
        safe_list,
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function solve_abstract_problem(abstract_problem::PR.OptimalControlProblem)
    abstract_controller = CO.NewControllerList()
    CO.compute_controller_reach!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.target_set,
    )
    return abstract_controller
end

function solve_abstract_problem(abstract_problem::PR.SafetyProblem)
    abstract_controller = CO.NewControllerList()
    CO.compute_controller_safe!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.safe_set,
    )
    return abstract_controller
end

function solve_concrete_problem(abstract_system, abstract_controller)
    function concrete_controller(x; param = false)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        upos = SY.get_upos_by_symbol(abstract_system, symbol)
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u
    end
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Build the abstraction
    abstract_system = build_abstraction(
        optimizer.concrete_problem.system,
        optimizer.state_grid,
        optimizer.input_grid,
    )
    optimizer.abstract_system = abstract_system
    # Build the abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem
    # Solve the abstract problem
    abstract_controller = solve_abstract_problem(abstract_problem)
    optimizer.abstract_controller = abstract_controller
    # Solve the concrete problem
    optimizer.concrete_controller =
        solve_concrete_problem(abstract_system, abstract_controller)

    optimizer.solve_time_sec = time() - t_ref
    return
end

end
