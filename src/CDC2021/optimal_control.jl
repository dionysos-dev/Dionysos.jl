module OptimalControl

using StaticArrays, Plots

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const U = Utils

using ..BranchAndBound
const BB = BranchAndBound

using ..Partition
const P = Partition

using ..AlternatingSimulation
const AS = AlternatingSimulation

using ..Lazy_abstraction
const LA = Lazy_abstraction

"""
-Udom::DomainList : the input domain
-_I_ : AB.HyperRectangle : necessary for the algo (must be as small as possible and contain x0)
-transition_cost(x,u) : function that calculates the transition cost from state x with input u
-functions : vector of problem-specific functions required
    compute_reachable_set(rect::AB.HyperRectangle,contsys,Udom::AB.DomainList)
    minimum_transition_cost(symmodel,contsys,source::Int,target::Int)
    post_image(symmodel,contsys,xpos,u)  : used in the A*
    pre_image(symmodel,contsys,xpos,u)   : used in the A*
-option : Vector of Bool to decide the branch and bound strategy (ex: indicate if you forbid the path to return to the previous big cell)
-ext : not used by default (could be used by the user)
"""
struct OptimalControlProblem{N,T,C,U,TC,CA,F<:Tuple{Vararg{Function}},E} <: BB.Abstract_BB_Problem
    # data of the problem
    x0::SVector{N,T}
    contsys::C
    periodic::Vector{Int}    # periodics dimensions
    periods::Vector{Float64} # periods
    T0::Vector{Float64}      # starting point of the period
    Udom::U
    transition_cost::TC
    # data of the algo
    q0::Int
    qT::Int
    coarse_abstraction::CA
    cells::Vector{P.Cell{N,T}}
    hx_medium::SVector{N,T}
    hx_fine::SVector{N,T}
    functions::F
    option::Vector{Bool}
    ext::E
end


function OptimalControlProblem(x0,_I_,_T_,contsys,periodic,periods,T0,Udom,transition_cost,partition,hx_medium,hx_fine,functions;option=[true], ext=nothing)
    q0,qT,coarse_abstraction,cells = P.Initialise(partition,contsys,Udom,_I_,_T_,functions[1],functions[2],periodic,periods,T0)
    return OptimalControlProblem(x0,contsys,periodic,periods,T0,Udom,transition_cost,q0,qT,coarse_abstraction,cells,hx_medium,hx_fine,functions,option,ext)
end

function BB.check_trivial_infeasibility(prob::OptimalControlProblem)
    return !(prob.cells[prob.q0].lower_bound < Inf)
end

function BB.get_first_instance(prob::OptimalControlProblem)
    return [prob.q0]
end

function build_medium_abstraction!(prob,cell)
    hx = prob.hx_medium
    Xdom = D.GeneralDomainList(hx;periodic=prob.periodic,periods=prob.periods,T0=prob.T0)
    AB.add_set!(Xdom, cell.hyperrectangle , AB.OUTER)
    for i in cell.outneighbors
        for rec in cell.local_target_set[i]
            AB.add_set!(Xdom,rec, AB.OUTER)
        end
    end
    symmodel = AB.NewSymbolicModelListList(Xdom, prob.Udom)
    problem = AS.symmodelProblem(symmodel,prob.contsys,prob.functions[1],prob.functions[2],AS.get_possible_transitions_2)
    autom = AS.build_alternating_simulation(problem)
    symmodel = AB.with_automaton(symmodel, autom)
    cell.medium_abstraction = symmodel
end
function build_heuristic!(cell,from)
    heuristics = cell.heuristics
    symmodel = cell.medium_abstraction
    _I_L = cell.local_init_set[from]
    initlist = U.get_symbols(symmodel,_I_L,AB.OUTER)
    heuristic = AS.build_heuristic(symmodel,initlist)
    heuristics[from] = heuristic
    #fig = plot(aspect_ratio = 1,legend = false)
    #AS.plot_heuristic!(heuristic,opacity=0.2,color=:red)
    #display(fig)
end


function BB.compute_lower_bound!(prob::OptimalControlProblem, node::BB.Node)
    path = node.elem
    cells = prob.cells
    curr_cell = cells[path[end]]
    from = length(path) == 1 ? -2 : path[end-1]
    cost = 0.0
    if curr_cell.index == prob.qT
        if curr_cell.medium_abstraction==nothing
            build_medium_abstraction!(prob,curr_cell)
        end
        if !haskey(curr_cell.heuristics,from)
            build_heuristic!(curr_cell,from)
        end
    end
    if length(path) > 1
        prev_cell = cells[path[end-1]]
        from_from = length(path) == 2 ? -2 : path[end-2]
        #if the alternating simulation is not yet computed
        if prev_cell.medium_abstraction==nothing
            build_medium_abstraction!(prob,prev_cell)
        end
        # build the heuristic if necessary
        if !haskey(prev_cell.heuristics,from_from)
            build_heuristic!(prev_cell,from_from)
        end

        heuristic = prev_cell.heuristics[from_from]
        _T_L = prev_cell.local_target_set[path[end]]
        cost = node.parent.lower_bound + AS.get_min_value_heurisitic(heuristic,_T_L)
        if prev_cell.index != prob.qT
            cost -= prev_cell.lower_bound
        else
            heuristic = prev_cell.heuristics[from_from]
            _T_L = prev_cell.local_target_set[-1]
            cost -= AS.get_min_value_heurisitic(heuristic,_T_L)
        end
    end

    if curr_cell.index != prob.qT
        cost += curr_cell.lower_bound
    else
        heuristic = curr_cell.heuristics[from]
        _T_L = curr_cell.local_target_set[-1]
        cost += AS.get_min_value_heurisitic(heuristic,_T_L)
    end
    #println(a)
    println("cost:  ",path," -> ",cost)
    node.lower_bound = cost#+ 5 # where 1 is the cost of the transition
    #println("end a lower bound")
end

function build_fine_abstraction!(prob,cell)
    hx = prob.hx_fine
    Xdom = D.GeneralDomainList(D.build_grid_in_rec(cell.hyperrectangle,hx),periodic=prob.periodic,periods=prob.periods,T0=prob.T0)
    #Xdom = D.GeneralDomainList(hx;periodic=prob.periodic,periods=prob.periods)
    AB.add_set!(Xdom, cell.hyperrectangle , AB.INNER)
    for i in cell.outneighbors
        for rec in cell.local_target_set[i]
            AB.add_set!(Xdom,rec, AB.OUTER) #INNER
        end
    end
    #U.plot_domain!(Xdom,opacity=0.4,color=:red)
    symmodel = AB.NewSymbolicModelListList(Xdom, prob.Udom, Set{NTuple{3,Int}})
    cell.fine_abstraction = (symmodel,nothing)
end
function build_controller!(prob,cell,prev,next)
    println("controller")
    (symmodel, transitions_previously_added) = cell.fine_abstraction
    _I_L = cell.local_init_set[prev]
    initlist = U.get_symbols(symmodel,_I_L,AB.OUTER)
    _T_L = cell.local_target_set[next]
    targetlist = U.get_symbols(symmodel,_T_L,AB.INNER)
    println(prev)
    println(next)
    println("length:")
    println(length(initlist))
    println(length(targetlist))
    heuristic = cell.heuristics[prev]

    problem,success = LA.compute_controller(symmodel, prob.contsys, initlist, targetlist,prob.transition_cost, prob.functions[4], prob.functions[3], h2, heuristic_data=heuristic,transitions_previously_added=transitions_previously_added)
    cell.fine_abstraction = (symmodel,problem.transitions_previously_added)
    cell.controllers[(prev,next)] = problem.contr
    println("plot")
    fig = plot(aspect_ratio = 1,legend = false)
    LA.plot_result!(problem)
    display(fig)
end

function BB.compute_upper_bound!(prob::OptimalControlProblem, node::BB.Node)
    path = node.elem
    cells = prob.cells
    if path[end]!=prob.qT
        return
    end
    #Threads.@threads for i=1:length(path)  # can be improve with a better repartition of the load on the threads...
    #    index = path[i]
    for (i,index) in enumerate(path)
        cell = cells[index]
        # if the abstraction (or any part of it) is not yet build
        if cell.fine_abstraction == nothing
            build_fine_abstraction!(prob,cell)
        end

        prev = i==1 ? -2 : path[i-1]
        next = i==length(path) ? -1 : path[i+1]
        controllers = cells[index].controllers
        # build the controller if not already computed
        if !haskey(controllers,(prev,next))
            build_controller!(prob,cell,prev,next)
        end
    end
    # compute the upper bound
    (traj, cost, succeed) = simulate_trajectory(prob, path)
    node.upper_bound = succeed ? cost : -Inf
    node.sol = path
    fig = plot(aspect_ratio = 1,legend = false)
    U.plot_domain!(prob.coarse_abstraction.Xdom,opacity=0.15,color=:blue)
    print_trajectory!(traj)
    display(fig)
end

function BB.expand(prob::OptimalControlProblem, node::BB.Node)
    path = node.elem
    current_cell_index = path[end]
    current_cell = prob.cells[current_cell_index]

    children = []
    for neighbor_cell in current_cell.outneighbors
        new_path = [path...,neighbor_cell]
        if prob.option[1] && (f1(reverse(new_path)) || (length(path) >= 2 && neighbor_cell == path[end-1]) || f2(new_path))
            continue
        end
        push!(children,BB.Node(new_path,node.depth+1,parent=node,lower_bound=-Inf, upper_bound=Inf,ext=prob)) #prob : temporary (to change)
    end
    return children
end

function simulate_trajectory(prob::OptimalControlProblem, path)
    x0 = prob.x0
    traj = []
    #push!(traj,x0)
    cost = 0.0
    for (i,index) in enumerate(path)
        cell = prob.cells[index]
        (symmodel,tab) = cell.fine_abstraction
        prev = i==1 ? -2 : path[i-1]
        next = i==length(path) ? -1 : path[i+1]
        contr = cell.controllers[(prev,next)]

        _T_L = cell.local_target_set[next]
        targetlist = U.get_symbols(symmodel,_T_L,AB.INNER)
        (local_traj,local_cost,succeed) = simulate_local_trajectory(prob.contsys, symmodel, contr, x0, targetlist; transition_cost=prob.transition_cost)
        #traj = [traj..., local_traj[2:end]...]
        traj = [traj..., local_traj[1:end-1]...]
        cost += local_cost
        if !succeed
            return (traj,cost,false)
        end
        x0,u = local_traj[end]
    end
    return (traj,cost,true)
end
function simulate_local_trajectory(contsys, symmodel, contr, x0, targetlist; transition_cost)
    traj = []
    cost = 0.0
    while true
        #push!(traj,x0)
        xpos = AB.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return (traj,cost,false)
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        if source ∈ targetlist
            push!(traj,(x0,nothing))
            break
        end
        symbollist = AB.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return (traj,cost,false)
        end
        symbol = first(symbollist)[1]

        upos = AB.get_upos_by_symbol(symmodel, symbol)
        u = AB.get_coord_by_pos(symmodel.Udom.grid, upos)
        cost += transition_cost(x0,u)

        push!(traj,(x0,u))
        x0 = contsys.sys_map(x0, u, contsys.tstep)
        x0 = D.set_in_period_coord(symmodel.Xdom,x0)
    end
    return (traj,cost,true)
end

function print_trajectory!(traj)
    for i=1:length(traj)-1
        x1,u1 = traj[i]
        x2,u2 = traj[i+1]
        plot!([x1[1],x2[1]], [x1[2],x2[2]],color =:red,linewidth = 2)
        if i>1
            scatter!([x1[1]],[x1[2]],color =:red,markersize=2)
        end
    end
    x1,u1 = traj[1]
    x2,u2 = traj[end]
    scatter!([x1[1]],[x1[2]],color =:green,markersize=3)
    scatter!([x2[1]],[x2[2]],color =:yellow,markersize=3)
end

# heurisic: Manhattan distance (could be discarded later)
function h1(node::LA.S.Node,problem::LA.LazyAbstraction)
    source = node.state.source
    tar = problem.goal[1].source
    symmodel = problem.symmodel
    xpos = AB.get_xpos_by_state(symmodel, source)
    xtar = AB.get_xpos_by_state(symmodel, tar)
    return (abs(xpos[1]-xtar[1]) + abs(xpos[2]-xtar[2]))
end

function h2(node::LA.S.Node,problem::LA.LazyAbstraction)
    source = node.state.source
    symmodel = problem.symmodel
    xpos = AB.get_xpos_by_state(symmodel, source)
    x = AB.get_coord_by_pos(symmodel.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = AB.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = AB.get_state_by_xpos(symmodel2, xpos2)
    return heuristic.dists[source2]
end

## heuristic function to eliminate nodes without considering them in the expand function

# return true if Vector A starts with two identical consecutive sequences
# used to impose that a trajectory can not make a cycle in the large cell graph
function f1(A)
    n = length(A)
    i = 1
    j = 2
    while j <= n && n-j>=j-2
        if A[j] == A[1]
            if j == 2
                return true
            end
            for k=1:j-2
                if A[i+k] != A[j+k]
                    break
                end
                if k == j-2
                    return true
                end
            end
        end
        j+=1
    end
    return false
end

# check if the last element is already in the array,
# used to impose that a trajectory can not visits twice the same large cell
function f2(A)
    last = A[end]
    for i=1:length(A)-1
        if(A[i]==last)
            return true
        end
    end
    return false
end
end #end module
