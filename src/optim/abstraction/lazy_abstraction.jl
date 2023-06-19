export LazyAbstraction

module LazyAbstraction
using JuMP
using LinearAlgebra, IntervalArithmetic, Random, Plots
Random.seed!(0)

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

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
    concrete_bell_fun::Union{Nothing, Any}

    lazy_search_problem::Union{Nothing, Any}

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
        )
    end
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "concrete_problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function set_abstract_system(concrete_system, hx, Ugrid)
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

    Udom = DO.DomainList(Ugrid)
    DO.add_set!(Udom, concrete_system.U, DO.OUTER)

    abstract_system = SY.LazySymbolicModel(Xdom, Udom)
    return abstract_system
end

function build_abstract_problem(concrete_problem::PR.OptimalControlProblem, abstract_system)
    initlist = SY.get_symbol(abstract_system, concrete_problem.initial_set, DO.OUTER)
    targetlist = SY.get_symbol(abstract_system, concrete_problem.target_set, DO.INNER)
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
function build_heuristic_data(concrete_problem, abstract_system, compute_reachable_set)
    function get_transitions(symmodel, sys, source)
        return SY.get_transitions_1(symmodel, sys, source, compute_reachable_set)
    end

    function minimum_transition_cost(symmodel, contsys, source, target)
        return 1.0
    end
    hx = [1.0, 1.0] * 1.5
    concrete_system = concrete_problem.system
    X = concrete_system.X.A
    Xdom = DO.GeneralDomainList(
        hx;
        periodic = concrete_system.periodic,
        periods = concrete_system.periods,
        T0 = concrete_system.T0,
    )
    DO.add_set!(Xdom, X, DO.OUTER)
    symmodel = SY.symmodelAS(
        Xdom,
        abstract_system.Udom,
        concrete_system,
        minimum_transition_cost,
        get_transitions,
    )
    initlist = SY.get_symbol(symmodel, concrete_problem.initial_set, DO.OUTER)
    heuristic_data = SY.build_heuristic(symmodel, initlist)
    return heuristic_data
end

function h1(node::UT.Node, problem)
    source = node.state.source
    abstract_system = problem.abstract_system
    xpos = SY.get_xpos_by_state(abstract_system, source)
    x = DO.get_coord_by_pos(abstract_system.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = DO.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = SY.get_state_by_xpos(symmodel2, xpos2)[1]
    return heuristic.dists[source2]
end

function set_optimizer!(
    optimizer::Optimizer,
    concrete_problem,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    hx,
    Ugrid;
    transitions_previously_added = nothing,
)
    optimizer.concrete_problem = concrete_problem
    concrete_system = concrete_problem.system
    optimizer.concrete_system = concrete_system

    # Initialize the abstract system
    abstract_system = set_abstract_system(concrete_system, hx, Ugrid)
    optimizer.abstract_system = abstract_system

    # Build the abstract specifications
    abstract_problem = build_abstract_problem(concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    # Build the abstraction-based heuristic
    heuristic_data =
        build_heuristic_data(concrete_problem, abstract_system, compute_reachable_set)

    # Build the search problem
    lazy_search_problem = LazySearchProblem(
        abstract_problem,
        concrete_problem,
        pre_image,
        post_image,
        compute_reachable_set,
        h1,
        heuristic_data,
        transitions_previously_added,
        maxIter,
    )
    optimizer.lazy_search_problem = lazy_search_problem

    return
end

function build_abstraction(optimizer::Optimizer)
    lazy_search_problem = optimizer.lazy_search_problem
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
        source = SY.get_state_by_xpos(abstract_system, xpos)[1]
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return
        end
        symbol = first(symbollist)[1]
        upos = SY.get_upos_by_symbol(abstract_system, symbol)
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u
    end
    return concrete_controller
end

function MOI.optimize!(optimizer::Optimizer)
    # Co-design the abstract system and the abstract controller
    build_abstraction(optimizer)
    lazy_search_problem = optimizer.lazy_search_problem
    abstract_system = optimizer.abstract_system
    abstract_controller = lazy_search_problem.contr
    optimizer.abstract_controller = abstract_controller
    optimizer.concrete_controller =
        get_concrete_controller(abstract_system, abstract_controller)

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
    contr = CO.NewControllerList()
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
    upos = SY.get_upos_by_symbol(problem.abstract_system, action)
    u = DO.get_coord_by_pos(problem.abstract_system.Udom.grid, upos)

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
        problem.costs[i] = 0.0
        problem.controllable[i] = false
    end
end

function update_abstraction!(successors, problem::LazySearchProblem, source)
    abstract_system = problem.abstract_system
    concrete_system = problem.concrete_system
    Xdom = abstract_system.Xdom
    Udom = abstract_system.Udom
    nsym = abstract_system.autom.nsymbols

    xpos = SY.get_xpos_by_state(abstract_system, source)
    x = DO.get_coord_by_pos(Xdom.grid, xpos)
    for upos in DO.enum_pos(Udom)
        symbol = SY.get_symbol_by_upos(abstract_system, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)

        ns1 = abstract_system.autom.nstates
        cells = problem.pre_image(abstract_system, concrete_system, xpos, u)
        _update_cache!(problem, ns1, abstract_system.autom.nstates, nsym)
        for cell in cells
            if !problem.controllable[cell]
                if !problem.transitions_added[cell, symbol]
                    # add transitions for input u starting from cell if it has not be done yet
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

## printing
function rectangle(c, r)
    return Shape(
        c[1] .- r[1] .+ [0, 2 * r[1], 2 * r[1], 0],
        c[2] .- r[2] .+ [0, 0, 2 * r[2], 2 * r[2]],
    )
end

function plot_result!(problem; dims = [1, 2], x0 = nothing)
    println()
    println("Plotting")
    targetlist = [init.source for init in problem.initial]
    initlist = [goal.source for goal in problem.goal]
    contsys = problem.concrete_system
    contr = problem.contr
    abstract_system = problem.abstract_system
    domain = abstract_system.Xdom
    grid = domain.grid
    h = grid.h[dims]

    # states for which transisitons have been computed for at least one input
    dict = Dict{NTuple{2, Int}, Any}()
    for k in 1:(abstract_system.autom.nstates)
        if any(u -> problem.transitions_added[k, u], 1:(problem.transitions_added.num_rows))
            pos = SY.get_xpos_by_state(abstract_system, k)
            if !haskey(dict, pos[dims])
                dict[pos[dims]] = true
                center = DO.get_coord_by_pos(grid, pos)
                Plots.plot!(rectangle(center[dims], h ./ 2); opacity = 0.2, color = :yellow)
            end
        end
    end

    # controllable state
    dict = Dict{NTuple{2, Int}, Any}()
    for (cell, symbol) in contr.data
        pos = SY.get_xpos_by_state(abstract_system, cell)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims], h ./ 2); opacity = 0.3, color = :blue)
        end
    end

    # states selected by A* to compute their pre-image
    dict = Dict{NTuple{2, Int}, Any}()
    for state in Base.keys(problem.closed)
        pos = SY.get_xpos_by_state(abstract_system, state.source)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims], h ./ 2); opacity = 0.5, color = :blue)
        end
    end

    # initial set
    dict = Dict{NTuple{2, Int}, Any}()
    for s in initlist
        pos = SY.get_xpos_by_state(abstract_system, s)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims], h ./ 2); opacity = 0.4, color = :green)
        end
    end

    # target set
    dict = Dict{NTuple{2, Int}, Any}()
    for s in targetlist
        pos = SY.get_xpos_by_state(abstract_system, s)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims], h ./ 2); opacity = 0.5, color = :red)
        end
    end

    # plot a trajectory
    if x0 != nothing
        (traj, success) = trajectory_reach(contsys, abstract_system, contr, x0, targetlist)
        plot_trajectory!(abstract_system, traj; dims = dims)
    end
end

function trajectory_reach(
    contsys,
    abstract_system,
    contr,
    x0,
    targetlist;
    randchoose = false,
)
    traj = []
    while true
        x0 = DO.set_in_period_coord(abstract_system.Xdom, x0)
        push!(traj, x0)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x0)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("Trajectory out of domain")
            return (traj, false)
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)[1]
        if source ∈ targetlist
            break
        end
        symbollist = UT.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return (traj, false)
        end
        if randchoose
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end

        upos = SY.get_upos_by_symbol(abstract_system, symbol)
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return (traj, true)
end

function plot_trajectory!(abstract_system, traj; dims = [1, 2])
    domain = abstract_system.Xdom
    grid = domain.grid
    k = dims[1]
    l = dims[2]
    for i in 1:(length(traj) - 1)
        Plots.plot!(
            [traj[i][k], traj[i + 1][k]],
            [traj[i][l], traj[i + 1][l]];
            color = :red,
            linewidth = 2,
        )
        if i > 1
            Plots.scatter!([traj[i][k]], [traj[i][l]]; color = :red, markersize = 2)
        end
    end
    Plots.scatter!([traj[1][k]], [traj[1][l]]; color = :green, markersize = 3)
    return Plots.scatter!([traj[end][k]], [traj[end][l]]; color = :yellow, markersize = 3)
end

end
