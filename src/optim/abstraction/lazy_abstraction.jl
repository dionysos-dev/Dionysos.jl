export LazyAbstraction

module LazyAbstraction
import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
using LinearAlgebra, JuMP, IntervalArithmetic, Random

Random.seed!(0)

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

"""
    Optimizer{T} <: MOI.AbstractOptimizer

Abstraction-based solver for which the hyper-rectangular abstraction and the controller are co-designed to reduce the computation cost of the abstraction.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_system::Union{Nothing, Any}
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, Any}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Union{Nothing, Any}
    abstract_lyap_fun::Union{Nothing, Any}
    concrete_lyap_fun::Union{Nothing, Any}
    abstract_bell_fun::Union{Nothing, Any}
    concrete_bell_fun::Union{Nothing, Any} # from the initial set
    abstract_system_heuristic::Union{Nothing, Any}
    lyap_fun::Union{Nothing, Any}
    bell_fun::Union{Nothing, Any}
    heuristic_data::Union{Nothing, Any}

    lazy_search_problem::Union{Nothing, Any}
    solved::Union{Nothing, Bool}
    param::Union{Nothing, Any}
    solve_time_sec::T

    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            nothing,
            0.0,
        )
    end
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "concrete_problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function set_abstract_system(concrete_system, hx, Udom; Xdom = nothing)
    if Xdom === nothing
        X = concrete_system.X.A
        obstacles = concrete_system.X.B.sets
        d = DO.RectangularObstacles(X, obstacles)
        Xdom = DO.GeneralDomainList(
            hx;
            elems = d,
            periodic = concrete_system.periodic,
            periods = concrete_system.periods,
            T0 = concrete_system.T0,
        )
    end

    # Udom = DO.DomainList(Ugrid)
    # DO.add_set!(Udom, concrete_system.U, DO.OUTER)
    abstract_system = SY.LazySymbolicModelList(Xdom, Udom)
    return abstract_system
end

function build_abstract_problem(concrete_problem::PR.OptimalControlProblem, abstract_system)
    initlist =
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, DO.OUTER)
    targetlist =
        SY.get_states_from_set(abstract_system, concrete_problem.target_set, DO.INNER)
    initial_set = [State(init) for init in initlist]
    target_set = [State(tar) for tar in targetlist]
    return PR.OptimalControlProblem(
        abstract_system,
        initial_set,
        target_set,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

# Construct the heuristic to guide the A* exploration
function build_abstract_system_heuristic(
    concrete_system,
    Udom,
    compute_reachable_set,
    minimum_transition_cost,
    hx_heuristic;
    Xdom = nothing,
)
    function get_transitions(symmodel, sys, source)
        return SY.get_transitions_1(symmodel, sys, source, compute_reachable_set)
    end
    if Xdom === nothing
        X = concrete_system.X.A
        Xdom = DO.GeneralDomainList(
            hx_heuristic;
            periodic = concrete_system.periodic,
            periods = concrete_system.periods,
            T0 = concrete_system.T0,
        )
        DO.add_set!(Xdom, X, DO.OUTER)
    end

    symmodel =
        SY.symmodelAS(Xdom, Udom, concrete_system, minimum_transition_cost, get_transitions)
    return symmodel
end

function build_heuristic(abstract_system_heuristic, initial_set)
    initlist = SY.get_states_from_set(abstract_system_heuristic, initial_set, DO.OUTER)
    heuristic_data = SY.build_heuristic(abstract_system_heuristic, initlist)
    return heuristic_data
end

# Default abstract heuristic
function h_abstract(node::UT.Node, problem)
    source = node.state.source
    abstract_system = problem.abstract_system
    xpos = SY.get_xpos_by_state(abstract_system, source)
    x = DO.get_coord_by_pos(abstract_system.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = DO.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = SY.get_state_by_xpos(symmodel2, xpos2)
    return heuristic.dists[source2]
end

function set_optimizer_parameters!(
    optimizer::Optimizer,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    minimum_transition_cost,
    hx,
    hx_heuristic;
    γ = 1.0,
    h_user = (node, prob) -> 0.0,
    transitions_previously_added = nothing,
)
    optimizer.param = Dict()
    optimizer.param[:maxIter] = maxIter
    optimizer.param[:pre_image] = pre_image
    optimizer.param[:post_image] = post_image
    optimizer.param[:compute_reachable_set] = compute_reachable_set
    optimizer.param[:minimum_transition_cost] = minimum_transition_cost
    optimizer.param[:hx] = hx
    optimizer.param[:hx_heuristic] = hx_heuristic
    optimizer.param[:γ] = γ
    optimizer.param[:h_user] = h_user
    return optimizer.param[:transitions_previously_added] = transitions_previously_added
end

function set_optimizer!(
    optimizer::Optimizer,
    concrete_problem,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    minimum_transition_cost,
    hx_heuristic,
    hx,
    Udom;
    abstract_system_heuristic = nothing,
    γ = 1.0,
    h_user = (node, prob) -> 0.0,
    transitions_previously_added = nothing,
    Xdom = nothing,
)
    set_optimizer_parameters!(
        optimizer,
        maxIter,
        pre_image,
        post_image,
        compute_reachable_set,
        minimum_transition_cost,
        hx,
        hx_heuristic;
        γ = γ,
        h_user = h_user,
        transitions_previously_added = transitions_previously_added,
    )

    optimizer.concrete_problem = concrete_problem
    concrete_system = concrete_problem.system
    optimizer.concrete_system = concrete_system

    # Initialize the abstract system
    abstract_system = set_abstract_system(concrete_system, hx, Udom; Xdom = Xdom)
    optimizer.abstract_system = abstract_system

    # Build the abstract specifications
    abstract_problem = build_abstract_problem(concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    # Set the abstract system for heuristic if given
    optimizer.abstract_system_heuristic = abstract_system_heuristic

    return
end

function set_optimizer!(
    optimizer::Optimizer,
    concrete_problem,
    Ugrid;
    abstract_system_heuristic = nothing,
    Xdom = nothing,
)
    return set_optimizer!(
        optimizer,
        concrete_problem,
        optimizer.param[:maxIter],
        optimizer.param[:pre_image],
        optimizer.param[:post_image],
        optimizer.param[:compute_reachable_set],
        optimizer.param[:minimum_transition_cost],
        optimizer.param[:hx_heuristic],
        optimizer.param[:hx],
        Ugrid;
        γ = optimizer.param[:γ],
        h_user = optimizer.param[:h_user],
        abstract_system_heuristic = abstract_system_heuristic,
        transitions_previously_added = nothing,
        Xdom = Xdom,
    )
end

function Base.copy(optimizer::Optimizer)
    new_optimizer = Optimizer()
    new_optimizer.param = copy(optimizer.param)
    return new_optimizer
end

function build_abstraction(optimizer::Optimizer)
    # Define the heuristic
    γ = optimizer.param[:γ]
    h_user = optimizer.param[:h_user]
    h(node::UT.Node, problem) = γ * max(h_abstract(node, problem), h_user(node, problem))

    # Build the search problem
    lazy_search_problem = LazySearchProblem(
        optimizer.abstract_problem,
        optimizer.concrete_problem,
        optimizer.param[:pre_image],
        optimizer.param[:post_image],
        optimizer.param[:compute_reachable_set],
        h,
        optimizer.heuristic_data,
        optimizer.param[:transitions_previously_added],
        optimizer.param[:maxIter],
    )
    optimizer.lazy_search_problem = lazy_search_problem

    node, nb = UT.astar_graph_search(lazy_search_problem, lazy_search_problem.h)
    println(
        "\nnumber of transitions created: ",
        length(lazy_search_problem.abstract_system.autom.transitions),
    )
    if node === nothing
        println("compute_controller_reach! terminated without covering init set")
        return optimizer, false
    end
    println("compute_controller_reach! terminated with success")
    return
end

function get_concrete_controller(abstract_system, abstract_controller)
    function concrete_controller(x)
        x = DO.set_in_period_coord(abstract_system.Xdom, x)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("Trajectory out of domain")
            return
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return
        end
        symbol = first(symbollist)[1]
        u = SY.get_concrete_input(abstract_system, symbol)
        return u
    end
    return concrete_controller
end

function build_abstract_lyap_fun(lyap)
    return abstract_lyap_fun(state) = lyap[state]
end

function build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)
    function concrete_lyap_fun(x)
        state = SY.get_abstract_state(abstract_system, x)
        return abstract_lyap_fun(state)
    end
    return concrete_lyap_fun
end

function build_abstract_bell_fun(bell)
    return abstract_bell_fun(state) = bell[state]
end

function build_concrete_bell_fun(abstract_system_heuristic, abstract_bell_fun)
    function concrete_bell_fun(x)
        state = SY.get_abstract_state(abstract_system_heuristic, x)
        return abstract_bell_fun(state)
    end
    return concrete_bell_fun
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    concrete_problem = optimizer.concrete_problem
    # Build the abstraction-based heuristic
    if optimizer.abstract_system_heuristic === nothing
        abstract_system_heuristic = build_abstract_system_heuristic(
            concrete_problem.system,
            optimizer.abstract_system.Udom,
            optimizer.param[:compute_reachable_set],
            optimizer.param[:minimum_transition_cost],
            optimizer.param[:hx_heuristic],
        )
        optimizer.abstract_system_heuristic = abstract_system_heuristic
    end
    # Build the Bellman-like value funtion
    heuristic_data =
        build_heuristic(optimizer.abstract_system_heuristic, concrete_problem.initial_set)
    optimizer.heuristic_data = heuristic_data
    optimizer.bell_fun =
        Dict(state => bell for (state, bell) in enumerate(heuristic_data.dists))

    # Co-design the abstract system and the abstract controller
    build_abstraction(optimizer)
    lazy_search_problem = optimizer.lazy_search_problem
    abstract_system = optimizer.abstract_system
    abstract_controller = lazy_search_problem.contr
    optimizer.abstract_controller = abstract_controller
    optimizer.concrete_controller =
        get_concrete_controller(abstract_system, abstract_controller)
    # Construct Lyapunov-like function
    lyap_fun = Dict(state => lyap for (state, lyap) in enumerate(lazy_search_problem.costs))
    optimizer.lyap_fun = lyap_fun
    abstract_lyap_fun = build_abstract_lyap_fun(lyap_fun)
    optimizer.abstract_lyap_fun = abstract_lyap_fun
    optimizer.concrete_lyap_fun =
        build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)

    # Construct Bellman-like value function 
    abstract_bell_fun = build_abstract_bell_fun(optimizer.bell_fun)
    optimizer.abstract_bell_fun = abstract_bell_fun
    optimizer.concrete_bell_fun =
        build_concrete_bell_fun(optimizer.abstract_system_heuristic, abstract_bell_fun)

    optimizer.solved = true

    optimizer.solve_time_sec = time() - t_ref
    return
end

struct State
    source::Int
end
Base.:(==)(s1::State, s2::State) = s1.source == s2.source

mutable struct MutableMatrix{T, VT <: AbstractVector{T}}
    default::T
    num_rows::Int
    num_cols::Int
    data::VT
end
_fill(default::Bool, N) = default ? trues(N) : falses(N)
_fill(default, N) = fill(default, N)
function MutableMatrix(default, num_rows::Int, num_cols::Int)
    N = num_rows * num_cols
    data = _fill(default, N)
    return MutableMatrix(default, num_rows, num_cols, data)
end
function add_columns!(m::MutableMatrix, n)
    start = length(m.data) + 1
    m.num_cols += n
    resize!(m.data, m.num_cols * m.num_rows)
    for i in start:length(m.data)
        m.data[i] = m.default
    end
end
Base.getindex(m::MutableMatrix, col::Int, row::Int) = m.data[(col - 1) * m.num_rows + row]
function Base.setindex!(m::MutableMatrix{T}, val::T, col::Int, row::Int) where {T}
    return m.data[(col - 1) * m.num_rows + row] = val
end

mutable struct LazySearchProblem{T} <: UT.SearchProblem{T}
    initial::Vector{T}
    goal::Vector{T}
    abstract_system::Union{Nothing, Any}
    concrete_system::Union{Nothing, Any}
    pre_image::Union{Nothing, Any}
    post_image::Union{Nothing, Any}
    compute_reachable_set::Union{Nothing, Any}
    transition_cost::Union{Nothing, Any}
    h::Union{Nothing, Any}
    heuristic_data::Union{Nothing, Any} # extension for potential additionnal data for the heuristic function
    contr::Union{Nothing, Any}

    transitions_added::Union{Nothing, MutableMatrix{Bool, BitVector}} # could be an array or a dictionnary (to be added)
    num_targets_unreachable::Union{Nothing, MutableMatrix{Int, Vector{Int}}} # could be an array or a dictionnary (to be added)
    controllable::Union{Nothing, BitVector} # could be an array or a dictionnary (to be added)
    num_init_unreachable::Union{Nothing, Int}  # counter of the remaining non controllable init cells
    closed::Union{Nothing, Dict{T, Bool}} # only usefull for the printing (could be discard later)
    costs_temp::Union{Nothing, MutableMatrix{Float64, Vector{Float64}}} # array containing the current worse cost to reach the target, if the next input applied is symbol
    costs::Union{Nothing, Vector{Float64}} # vector containing the (worst) cost to reach the target set for each cell (necessary because of the pseudo non determinism) = Lyapunov function
    transitions_previously_added::Union{Nothing, MutableMatrix{Int, Vector{Int}}} # only necessary, if we need to reuse a partially computed symmodel
    # this array should contain the number of outgoing neighbors for
    # previously computed couple (cell,input) and -1 for the others.
    maxIter::Union{Nothing, Int}
end

function LazySearchProblem(
    abstract_problem,
    concrete_problem,
    pre_image,
    post_image,
    compute_reachable_set,
    h,
    heuristic_data,
    transitions_previously_added,
    maxIter,
)
    initial = abstract_problem.target_set
    goal = abstract_problem.initial_set
    abstract_system = abstract_problem.system
    concrete_system = concrete_problem.system

    transition_cost(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
    transitions_added =
        MutableMatrix(false, abstract_system.autom.nsymbols, abstract_system.autom.nstates)
    num_targets_unreachable =
        MutableMatrix(0, abstract_system.autom.nsymbols, abstract_system.autom.nstates)
    controllable = falses(abstract_system.autom.nstates)
    for init in initial
        controllable[init.source] = true
    end
    num_init_unreachable = length(goal)
    costs_temp =
        MutableMatrix(0.0, abstract_system.autom.nsymbols, abstract_system.autom.nstates)
    costs = zeros(Float64, abstract_system.autom.nstates)
    if transitions_previously_added === nothing
        transitions_previously_added =
            MutableMatrix(-1, abstract_system.autom.nsymbols, abstract_system.autom.nstates)
    end

    closed = nothing
    contr = UT.SortedTupleSet{2, NTuple{2, Int}}()
    return LazySearchProblem(
        initial,
        goal,
        abstract_system,
        concrete_system,
        pre_image,
        post_image,
        compute_reachable_set,
        transition_cost,
        h,
        heuristic_data,
        contr,
        transitions_added,
        num_targets_unreachable,
        controllable,
        num_init_unreachable,
        closed,
        costs_temp,
        costs,
        transitions_previously_added,
        maxIter,
    )
end

function UT.goal_test(problem::LazySearchProblem, state::State)
    if state.source in [s.source for s in problem.goal]
        if iszero(problem.num_init_unreachable -= 1)
            return true
        end
    end
    return false
end

# for the moment, one action costs 1
function UT.path_cost(problem::LazySearchProblem, c, state1::State, action, state2::State)
    source = state2.source
    pos = SY.get_xpos_by_state(problem.abstract_system, source)
    x = DO.get_coord_by_pos(problem.abstract_system.Xdom.grid, pos)
    u = SY.get_concrete_input(problem.abstract_system, action)

    problem.costs[source] += problem.transition_cost(x, u)
    return problem.costs[source]
end

function transitions!(source, symbol, u, abstract_system, concrete_system, post_image)
    xpos = SY.get_xpos_by_state(abstract_system, source)
    over_approx = post_image(abstract_system, concrete_system, xpos, u)
    translist = [(cell, source, symbol) for cell in over_approx]
    SY.add_transitions!(abstract_system.autom, translist)
    return length(over_approx)
end

function _update_cache!(problem::LazySearchProblem, ns1, ns2, nsym)
    ns2 == ns1 && return
    Δ = ns2 - ns1
    add_columns!(problem.transitions_added, Δ)
    add_columns!(problem.num_targets_unreachable, Δ)
    add_columns!(problem.transitions_previously_added, Δ)
    add_columns!(problem.costs_temp, Δ)
    resize!(problem.controllable, ns2)
    resize!(problem.costs, ns2)
    for i in (ns1 + 1):ns2
        problem.costs[i] = Inf
        problem.controllable[i] = false
    end
end

function update_abstraction!(successors, problem::LazySearchProblem, source)
    abstract_system = problem.abstract_system
    concrete_system = problem.concrete_system
    Udom = abstract_system.Udom
    nsym = abstract_system.autom.nsymbols

    xpos = SY.get_xpos_by_state(abstract_system, source)
    for symbol in SY.enum_inputs(abstract_system)
        u = SY.get_concrete_input(abstract_system, symbol)
        ns1 = abstract_system.autom.nstates
        cells = problem.pre_image(abstract_system, concrete_system, xpos, u)
        _update_cache!(problem, ns1, abstract_system.autom.nstates, nsym)
        for cell in cells
            if !problem.controllable[cell]
                if !problem.transitions_added[cell, symbol]
                    n_trans = 0
                    if problem.transitions_previously_added[cell, symbol] != -1
                        n_trans = problem.transitions_previously_added[cell, symbol]
                    else
                        ns1 = abstract_system.autom.nstates
                        n_trans = transitions!(
                            cell,
                            symbol,
                            u,
                            abstract_system,
                            concrete_system,
                            problem.post_image,
                        )
                        _update_cache!(problem, ns1, abstract_system.autom.nstates, nsym)
                        problem.transitions_previously_added[cell, symbol] = n_trans
                    end
                    problem.num_targets_unreachable[cell, symbol] = n_trans
                    problem.transitions_added[cell, symbol] = true
                end
                # check if the cell is really in the pre-image
                if (source, cell, symbol) in abstract_system.autom.transitions
                    problem.costs_temp[cell, symbol] =
                        max(problem.costs_temp[cell, symbol], problem.costs[source])
                    if iszero(problem.num_targets_unreachable[cell, symbol] -= 1)
                        problem.costs[cell] = problem.costs_temp[cell, symbol]
                        problem.controllable[cell] = true
                        push!(successors, (symbol, State(cell)))
                        UT.push_new!(problem.contr, (cell, symbol))
                    end
                end
            end
        end
    end
end

function UT.successor(problem::LazySearchProblem, state::State)
    successors = []
    update_abstraction!(successors, problem, state.source)
    return successors
end

@recipe function f(problem::LazySearchProblem; dims = [1, 2])
    targetlist = [init.source for init in problem.initial]
    initlist = [goal.source for goal in problem.goal]
    contr = problem.contr
    abstract_system = problem.abstract_system
    domain = abstract_system.Xdom
    grid = domain.grid
    legend := false
    # states for which transisitons have been computed for at least one input
    dict = Dict{NTuple{2, Int}, Any}()
    for k in 1:(abstract_system.autom.nstates)
        if any(u -> problem.transitions_added[k, u], 1:(problem.transitions_added.num_rows))
            pos = SY.get_xpos_by_state(abstract_system, k)
            if !haskey(dict, pos[dims])
                dict[pos[dims]] = true
                @series begin
                    opacity := 0.2
                    color := :yellow
                    return grid, pos
                end
            end
        end
    end

    # controllable state
    dict = Dict{NTuple{2, Int}, Any}()
    for (cell, symbol) in contr.data
        pos = SY.get_xpos_by_state(abstract_system, cell)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            @series begin
                opacity := 0.3
                color := :blue
                return grid, pos
            end
        end
    end

    # states selected by A* to compute their pre-image
    dict = Dict{NTuple{2, Int}, Any}()
    for state in Base.keys(problem.closed)
        pos = SY.get_xpos_by_state(abstract_system, state.source)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            @series begin
                opacity := 0.5
                color := :blue
                return grid, pos
            end
        end
    end

    # initial set
    dict = Dict{NTuple{2, Int}, Any}()
    for s in initlist
        pos = SY.get_xpos_by_state(abstract_system, s)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            @series begin
                opacity := 0.4
                color := :green
                return grid, pos
            end
        end
    end

    # target set
    dict = Dict{NTuple{2, Int}, Any}()
    for s in targetlist
        pos = SY.get_xpos_by_state(abstract_system, s)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            @series begin
                opacity := 0.5
                color := :red
                return grid, pos
            end
        end
    end
end

end
