export HierarchicalAbstraction

module HierarchicalAbstraction
import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series

using LinearAlgebra, JuMP, Graphs, SimpleWeightedGraphs, IntervalArithmetic
using Base.Threads

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const BB = UT.BranchAndBound

using ..LazyAbstraction

"""
    Optimizer{T} <: MOI.AbstractOptimizer
    
Abstraction-based solver for which the domain is initially partioned into coarse hyper-rectangular cells, which are iteratively locally smartly refined with respect to the control task.
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

    reference_local_optimizer::Union{Nothing, Any}
    hierarchical_problem::Union{Nothing, Any}
    optimizer_BB::Union{Nothing, Any}
    max_iter::Union{Nothing, Int}
    max_time::Union{Nothing, Any}
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

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end
function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function build_abstract_system(
    concrete_system,
    hx,
    Ugrid,
    compute_reachable_set,
    minimum_transition_cost,
)
    function get_transitions(symmodel, sys, source)
        return SY.get_transitions_1(symmodel, sys, source, compute_reachable_set)
    end
    X = concrete_system.X.A
    obstacles = concrete_system.X.B.sets

    d = DO.RectangularObstacles(X, obstacles)
    Xdom = DO.GeneralDomainList(
        hx;
        elems = d,
        periodic = concrete_system.periodic,
        periods = concrete_system.periods,
        T0 = concrete_system.T0,
        fit = true,
    )

    Udom = DO.DomainList(Ugrid)
    DO.add_set!(Udom, concrete_system.U, DO.OUTER)

    symmodel =
        SY.symmodelAS(Xdom, Udom, concrete_system, minimum_transition_cost, get_transitions)
    sub_symmodels = Dict()
    abstract_system = SY.HierarchicalSymbolicSystem(symmodel, sub_symmodels)
    return abstract_system
end

function build_abstract_problem(concrete_problem::PR.OptimalControlProblem, abstract_system)
    initial_set = SY.get_symbol(abstract_system, concrete_problem.initial_set, DO.OUTER)
    target_set = SY.get_symbol(abstract_system, concrete_problem.target_set, DO.OUTER)
    return PR.OptimalControlProblem(
        abstract_system,
        initial_set,
        target_set,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function set_optimizer!(
    optimizer::Optimizer,
    concrete_problem,
    hx_global,
    Ugrid,
    compute_reachable_set,
    minimum_transition_cost,
    reference_local_optimizer,
    max_iter,
    max_time;
    option = false,
)
    # Set the algorithm parameters
    optimizer.param = Dict()
    optimizer.param[:compute_reachable_set] = compute_reachable_set
    optimizer.param[:minimum_transition_cost] = minimum_transition_cost
    optimizer.param[:option] = option
    optimizer.max_iter = max_iter
    optimizer.max_time = max_time

    optimizer.concrete_problem = concrete_problem
    concrete_system = concrete_problem.system
    optimizer.concrete_system = concrete_system
    optimizer.reference_local_optimizer = reference_local_optimizer

    abstract_system = build_abstract_system(
        concrete_system,
        hx_global,
        Ugrid,
        compute_reachable_set,
        minimum_transition_cost,
    )
    optimizer.abstract_system = abstract_system
    abstract_problem = build_abstract_problem(concrete_problem, abstract_system)
    return optimizer.abstract_problem = abstract_problem
end

function simulate_trajectory(optimizer::Optimizer, x0)
    best_path = optimizer.optimizer_BB.best_sol
    return simulate_trajectory(optimizer.hierarchical_problem, best_path, x0)
end

function get_closed_loop_trajectory(optimizer::Optimizer, x0)
    x_traj, u_traj, c_traj = simulate_trajectory(optimizer, x0)
    control_trajectory = ST.Control_trajectory(ST.Trajectory(x_traj), ST.Trajectory(u_traj))
    return ST.Cost_control_trajectory(control_trajectory, ST.Trajectory(c_traj))
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Co-design the abstract system and the abstract controller
    hierarchical_problem = HierarchicalProblem(
        optimizer.concrete_problem,
        optimizer.abstract_problem,
        optimizer.param[:compute_reachable_set],
        optimizer.reference_local_optimizer;
        option = [optimizer.param[:option]],
    )
    optimizer.hierarchical_problem = hierarchical_problem
    optimizer_BB = BB.Optimizer(
        hierarchical_problem,
        optimizer.max_iter,
        optimizer.max_time;
        log_level = 2,
    )
    optimizer.optimizer_BB = optimizer_BB
    MOI.optimize!(optimizer_BB)
    optimizer.solved = BB.has_solution(optimizer_BB)

    optimizer.solve_time_sec = time() - t_ref

    return
end

#########################################################################################################################

mutable struct Cell
    state::Int
    outneighbors::Any
    inneighbors::Any
    local_target_set::Any
    local_init_set::Any
    heuristic_abstraction::Any
    heuristics::Any
    optimizers::Any
    upper_bound::Any
    lower_bound::Any
    rec::Any
    reachable_set::Any
end

mutable struct HierarchicalProblem <: BB.Abstract_BB_Problem
    concrete_problem::Any
    concrete_system::Any
    abstract_problem::Any # global abstract problem
    abstract_system::Any # global abstract system coarse_abstraction
    reference_local_optimizer::Any

    cells::Any
    transition_cost::Any
    x0::Any
    q0::Any
    qT::Any
    option::Any
    ext::Any
end

function compute_reachable_sets(concrete_system, abstract_system, compute_reachable_set)
    reachable_set_sets = Dict()
    symmodel = abstract_system.symmodel
    for state in SY.enum_states(abstract_system)
        pos = SY.get_xpos_by_state(abstract_system, state)
        rec = DO.get_rec(DO.get_grid(symmodel.Xdom), pos)
        reachable_set = compute_reachable_set(rec, concrete_system, symmodel.Udom)
        reachable_set_sets[state] = reachable_set
    end
    return reachable_set_sets
end

function compute_local_sets!(cells, q0, _I_, qT, _T_, periodic, periods, T0)
    cells[q0].local_init_set[-2] = []
    cells[qT].local_target_set[-1] = []
    for (state, cell) in cells
        for i in cell.inneighbors
            cell.local_init_set[i] = []
        end
        for i in cell.outneighbors
            cell.local_target_set[i] = []
        end
    end
    push!(cells[q0].local_init_set[-2], DO.set_rec_in_period(periodic, periods, T0, _I_)...)
    push!(
        cells[qT].local_target_set[-1],
        DO.set_rec_in_period(periodic, periods, T0, _T_)...,
    )
    for (state, cell) in cells
        RL = DO.set_rec_in_period(periodic, periods, T0, cell.reachable_set)
        for rec in RL
            for i in cell.outneighbors
                neigh = cells[i]
                I = Base.intersect(rec, neigh.rec)
                push!(neigh.local_init_set[cell.state], I)
                push!(cell.local_target_set[neigh.state], I)
            end
        end
    end
end

function build_cells(
    concrete_problem,
    abstract_problem,
    compute_reachable_set,
    q0,
    _I_,
    qT,
    _T_,
    reference_local_optimizer,
)
    concrete_system = concrete_problem.system
    abstract_system = abstract_problem.system
    reachable_set_sets =
        compute_reachable_sets(concrete_system, abstract_system, compute_reachable_set)
    heuristic = SY.build_heuristic(abstract_system.symmodel, [q0])
    symmodel = abstract_system.symmodel
    cells = Dict()
    for state in SY.enum_states(abstract_system)
        pos = SY.get_xpos_by_state(abstract_system, state)
        outneighbors = SimpleWeightedGraphs.outneighbors(symmodel.autom, state)
        inneighbors = SimpleWeightedGraphs.inneighbors(symmodel.autom, state)
        heuristic_abstraction = nothing
        heuristics = Dict{Int, Any}()
        lower_bound = heuristic.dists[state]
        upper_bound = Inf
        rec = DO.get_rec(symmodel.Xdom.grid, pos)
        reachable_set = reachable_set_sets[state]
        local_target_set = Dict{Int, Vector{typeof(_I_)}}()
        local_init_set = Dict{Int, Vector{typeof(_I_)}}()

        optimizers = Dict()
        cell = Cell(
            state,
            outneighbors,
            inneighbors,
            local_target_set,
            local_init_set,
            heuristic_abstraction,
            heuristics,
            optimizers,
            upper_bound,
            lower_bound,
            rec,
            reachable_set,
        )
        cells[state] = cell
    end
    compute_local_sets!(
        cells,
        q0,
        _I_,
        qT,
        _T_,
        concrete_system.periodic,
        concrete_system.periods,
        concrete_system.T0,
    )
    return cells
end

function HierarchicalProblem(
    concrete_problem,
    abstract_problem,
    compute_reachable_set,
    reference_local_optimizer;
    option = [],
    ext = [],
)
    concrete_system = concrete_problem.system
    abstract_system = abstract_problem.system
    q0 = abstract_problem.initial_set[1]
    qT = abstract_problem.target_set[1]
    _I_ = concrete_problem.initial_set
    _T_ = concrete_problem.target_set
    cells = build_cells(
        concrete_problem,
        abstract_problem,
        compute_reachable_set,
        q0,
        _I_,
        qT,
        _T_,
        reference_local_optimizer,
    )
    x0 = UT.get_center(concrete_problem.initial_set)
    return HierarchicalProblem(
        concrete_problem,
        concrete_system,
        abstract_problem,
        abstract_system,
        reference_local_optimizer,
        cells,
        abstract_problem.transition_cost,
        x0,
        q0,
        qT,
        option,
        ext,
    )
end

function BB.check_trivial_infeasibility(prob::HierarchicalProblem)
    return !(prob.cells[prob.q0].lower_bound < Inf)
end

function BB.get_first_instance(prob::HierarchicalProblem)
    return [prob.q0]
end

function build_heuristic_abstraction!(prob::HierarchicalProblem, cell::Cell)
    concrete_system = prob.concrete_system
    hx_heuristic = prob.reference_local_optimizer.param[:hx_heuristic]
    compute_reachable_set = prob.reference_local_optimizer.param[:compute_reachable_set]
    minimum_transition_cost = prob.reference_local_optimizer.param[:minimum_transition_cost]
    # Build Xdom
    Xdom = DO.GeneralDomainList(
        hx_heuristic;
        periodic = concrete_system.periodic,
        periods = concrete_system.periods,
        T0 = concrete_system.T0,
    )
    DO.add_set!(Xdom, cell.rec, DO.OUTER)
    for i in cell.outneighbors
        for rec in cell.local_target_set[i]
            DO.add_set!(Xdom, rec, DO.OUTER)
        end
    end
    v = 2.0
    rec = UT.HyperRectangle(
        cell.reachable_set.lb - SVector(v, v),
        cell.reachable_set.ub + SVector(v, v),
    )
    DO.add_set!(Xdom, rec, DO.OUTER)
    abstract_heuristic_system = LazyAbstraction.build_abstract_system_heuristic(
        concrete_system,
        prob.abstract_system.symmodel.Udom,
        compute_reachable_set,
        minimum_transition_cost,
        hx_heuristic;
        Xdom = Xdom,
    )
    cell.heuristic_abstraction = abstract_heuristic_system
    return cell.heuristic_abstraction
end

function build_heuristic!(cell::Cell, from)
    heuristic_abstraction = cell.heuristic_abstraction
    initial_set = cell.local_init_set[from]
    heuristic_data = LazyAbstraction.build_heuristic(heuristic_abstraction, initial_set[1])
    cell.heuristics[from] = heuristic_data
    return
end

function get_min_value_heurisitic(heuristic, subsetList)
    symmodel = heuristic.symmodel
    val = Inf
    for subset in subsetList
        posL = DO.get_subset_pos(symmodel.Xdom, subset, DO.OUTER)
        for pos in posL
            val = min(val, heuristic.dists[SY.get_state_by_xpos(symmodel, pos)])
        end
    end
    return val
end

function BB.compute_lower_bound!(prob::HierarchicalProblem, node::BB.Node)
    path = node.elem
    cells = prob.cells
    curr_cell = cells[path[end]]
    from = length(path) == 1 ? -2 : path[end - 1]
    cost = 0.0
    if curr_cell.state == prob.qT
        if curr_cell.heuristic_abstraction == nothing
            build_heuristic_abstraction!(prob, curr_cell)
        end
        if !haskey(curr_cell.heuristics, from)
            build_heuristic!(curr_cell, from)
        end
    end
    if length(path) > 1
        prev_cell = cells[path[end - 1]]
        from_from = length(path) == 2 ? -2 : path[end - 2]
        #if the alternating simulation is not yet computed
        if prev_cell.heuristic_abstraction == nothing
            build_heuristic_abstraction!(prob, prev_cell)
        end
        # build the heuristic if necessary
        if !haskey(prev_cell.heuristics, from_from)
            build_heuristic!(prev_cell, from_from)
        end

        heuristic = prev_cell.heuristics[from_from]
        _T_L = prev_cell.local_target_set[path[end]]
        cost = node.parent.lower_bound + get_min_value_heurisitic(heuristic, _T_L)
        if prev_cell.state != prob.qT
            cost -= prev_cell.lower_bound
        else
            heuristic = prev_cell.heuristics[from_from]
            _T_L = prev_cell.local_target_set[-1]
            cost -= get_min_value_heurisitic(heuristic, _T_L)
        end
    end

    if curr_cell.state != prob.qT
        cost += curr_cell.lower_bound
    else
        heuristic = curr_cell.heuristics[from]
        _T_L = curr_cell.local_target_set[-1]
        cost += get_min_value_heurisitic(heuristic, _T_L)
    end
    println("cost:  ", path, " -> ", cost)
    return node.lower_bound = cost
end

function set_local_optimizer!(
    cell,
    concrete_problem,
    abstract_system,
    reference_local_optimizer,
)
    optimizers = cell.optimizers
    concrete_system = concrete_problem.system
    v = 1.0
    rec = UT.HyperRectangle(
        cell.reachable_set.lb - SVector(v, v),
        cell.reachable_set.ub + SVector(v, v),
    )
    d = DO.RectangularObstacles(rec, [])
    Xdom = DO.GeneralDomainList(
        reference_local_optimizer.param[:hx];
        elems = d,
        periodic = concrete_system.periodic,
        periods = concrete_system.periods,
        T0 = concrete_system.T0,
        fit = true,
    )

    for (prev_state, _I_L) in cell.local_init_set
        for (next_state, _T_L) in cell.local_target_set
            local_concrete_problem = PR.OptimalControlProblem(
                concrete_problem.system,
                _I_L[1],
                _T_L[1],
                concrete_problem.state_cost,
                concrete_problem.transition_cost,
                concrete_problem.time,
            )
            local_optimizer = copy(reference_local_optimizer)
            LazyAbstraction.set_optimizer!(
                local_optimizer,
                local_concrete_problem,
                abstract_system.symmodel.Udom; #DO.get_grid(abstract_system.symmodel.Udom);
                Xdom = Xdom,
            )
            local_optimizer.abstract_system_heuristic = cell.heuristic_abstraction
            optimizers[(prev_state, next_state)] = local_optimizer
        end
    end
end

function BB.compute_upper_bound!(prob::HierarchicalProblem, node::BB.Node)
    global fig
    path = node.elem
    cells = prob.cells
    if path[end] != prob.qT
        return
    end
    n_threads = 4
    lk = ReentrantLock()
    Threads.nthreads() < n_threads
    Threads.@threads for i in 1:length(path)
        println("Iteration $i is running on thread $(Threads.threadid())")
        index = path[i]
        cell = cells[index]
        if isempty(cell.optimizers)
            set_local_optimizer!(
                cell,
                prob.concrete_problem,
                prob.abstract_system,
                prob.reference_local_optimizer,
            )
        end
        prev = i == 1 ? -2 : path[i - 1]
        next = i == length(path) ? -1 : path[i + 1]
        local_optimizer = cell.optimizers[(prev, next)]
        if !local_optimizer.solved
            MOI.optimize!(local_optimizer)
        end
    end
    (x_traj, u_traj, cost_traj) = simulate_trajectory(prob, path, prob.x0)
    succeed = true
    cost = sum(cost_traj)
    node.upper_bound = succeed ? cost : -Inf
    node.sol = path
    return
end

function BB.expand(prob::HierarchicalProblem, node::BB.Node)
    path = node.elem
    current_cell_index = path[end]
    current_cell = prob.cells[current_cell_index]

    childrens = []
    bool = true
    for neighbor_cell in current_cell.outneighbors
        new_path = [path..., neighbor_cell]
        if f1(new_path) || !bool
            continue
        end
        if prob.option[1]
            new_path = [path..., f3(length(path) + 1, prob)]
            bool = false
        end

        push!(
            childrens,
            BB.Node(
                new_path,
                node.depth + 1;
                parent = node,
                lower_bound = -Inf,
                upper_bound = Inf,
                ext = prob,
            ),
        )
    end
    return childrens
end

function simulate_trajectory(prob::HierarchicalProblem, path, x0)
    full_x_traj = [x0]
    full_u_traj = []
    full_cost_traj = []
    for (i, index) in enumerate(path)
        cell = prob.cells[index]
        prev = i == 1 ? -2 : path[i - 1]
        next = i == length(path) ? -1 : path[i + 1]
        local_optimizer = cell.optimizers[(prev, next)]
        concrete_problem = local_optimizer.concrete_problem
        concrete_system = concrete_problem.system
        concrete_controller = local_optimizer.concrete_controller
        reached(x) = x ∈ concrete_problem.target_set
        cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
        nstep = 100
        cost_control_trajectory = ST.get_closed_loop_trajectory(
            concrete_system.f_eval,
            concrete_controller,
            cost_eval,
            x0,
            nstep;
            stopping = reached,
            noise = false,
        )
        x_traj = cost_control_trajectory.control_trajectory.states.seq
        u_traj = cost_control_trajectory.control_trajectory.inputs.seq
        cost_traj = cost_control_trajectory.costs.seq
        append!(full_x_traj, x_traj[2:end])
        append!(full_u_traj, u_traj)
        append!(full_cost_traj, cost_traj)
        x0 = x_traj[end]
    end
    return (full_x_traj, full_u_traj, full_cost_traj)
end

function get_shortest_path(graph, start, target)
    results = Graphs.dijkstra_shortest_paths(graph, [target])
    parents = results.parents
    path = [start]
    i = start
    while parents[i] !== 0
        i = parents[i]
        push!(path, i)
    end
    return path
end

function get_shortest_path_abstract_system(abstract_system, abstract_problem)
    symmodel = abstract_system.symmodel
    autom = symmodel.autom
    start = abstract_problem.initial_set[1]
    target = abstract_problem.target_set[1]
    return path = get_shortest_path(autom, start, target)
end

function f1(tableau)
    occurences = Dict()
    for element in tableau
        if haskey(occurences, element)
            return true
        end
        occurences[element] = 1
    end
    return false
end

function f3(i, prob)
    new_path = [
        get_shortest_path_abstract_system(prob.abstract_system, prob.abstract_problem)...,
        28,
    ]
    new_path = [1, 2, 3, 4, 5, 11, 16, 21, 27, 28]
    new_path = [1, 2, 9, 15, 21, 27, 28]
    return new_path[i]
end

@recipe function f(prob::HierarchicalProblem; path = [], heuristic = false, fine = true)
    @series begin
        color := :blue
        return prob.abstract_system.symmodel.Xdom
    end
    if heuristic
        for (i, index) in enumerate(path)
            cell = prob.cells[index]
            prev = i == 1 ? -2 : path[i - 1]
            heuristic_data = cell.heuristics[prev]
            bell_fun =
                Dict(state => bell for (state, bell) in enumerate(heuristic_data.dists))
            @series begin
                arrowsB := false
                dims := [1, 2]
                cost := true
                lyap_fun := bell_fun
                return cell.heuristic_abstraction
            end
        end
    end
    if fine
        for (i, index) in enumerate(path)
            cell = prob.cells[index]
            prev = i == 1 ? -2 : path[i - 1]
            next = i == length(path) ? -1 : path[i + 1]
            local_optimizer = cell.optimizers[(prev, next)]
            @series begin
                return local_optimizer.lazy_search_problem
            end
        end
    end
end

end
