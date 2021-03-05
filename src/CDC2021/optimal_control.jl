module OptimalControl

using StaticArrays, Plots

using ..Abstraction
AB = Abstraction

using ..Utils
U = Utils

using ..BranchAndBound
BB = BranchAndBound

using ..Partition
P = Partition

using ..AlternatingSimulation
AS = AlternatingSimulation

using ..Lazy_abstraction
LA = Lazy_abstraction

"""
-Udom::DomainList : the input domain
-_I_ : AB.HyperRectangle : necessary for the algo (must be as small as possible and contain x0)
-transition_cost(x,u) : function that calculates the transition cost from state x with input u
-functions : vector of problem-specific functions needed
    compute_reachable_set(rect::AB.HyperRectangle,contsys,Udom::AB.DomainList)
    minimum_transition_cost(symmodel,contsys,source::Int,target::Int)
    post_image(symmodel,contsys,xpos,u)  : used in the A*
    pre_image(symmodel,contsys,xpos,u)   : used in the A*
-option : Bool to indicate if you forbid the path to return to the previous big cell
-ext : not used by default (could be used by the user)
"""
struct OptimalControlProblem <: BB.Abstract_BB_Problem
    q0::Int
    x0
    qT::Int
    _T_
    contsys
    Udom
    transition_cost
    cells
    option::Bool
    functions::Vector{Any}
    dim
    ext
end

function OptimalControlProblem(rectangles,x0,_I_,_T_,contsys,Udom,transition_cost,functions;option=true, ext=nothing)
    q0,qT,cells = P.Initialise(rectangles,contsys,Udom,_I_,_T_,functions[1],functions[2])
    dim = U.dims(rectangles[1])
    return OptimalControlProblem(q0,x0,qT,_T_,contsys,Udom,transition_cost,cells,option,functions,dim,ext)
end

function BB.check_trivial_infeasibility(prob::OptimalControlProblem)
    return !(prob.cells[prob.q0].lower_bound < Inf)
end

function BB.get_first_instance(prob::OptimalControlProblem)
    return [prob.q0]
end

function build_coarse_abstraction!(prob,cell)
    hx = SVector(1.0, 1.0)*0.5; x0 = @SVector zeros(prob.dim)
    Xgrid = AB.GridFree(x0,hx)
    Xdom = AB.DomainList(Xgrid)
    AB.add_set!(Xdom, cell.hyperrectangle , AB.OUTER)
    for i in cell.outneighbors
        AB.add_set!(Xdom,cell.local_target_set[i], AB.OUTER)
    end
    symmodel = AB.NewSymbolicModelListList(Xdom, prob.Udom)
    problem = AS.symmodelProblem(symmodel,prob.contsys,prob.functions[1],prob.functions[2],AS.get_possible_transitions_2)
    symmodel.autom = AS.build_alternating_simulation(problem)
    cell.coarse_abstraction = symmodel
end
function build_heuristic!(cell,from)
    heuristics = cell.heuristics
    symmodel = cell.coarse_abstraction
    _I_ = cell.local_init_set[from]
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    heuristic = AS.build_heuristic(symmodel,initlist)
    heuristics[from] = heuristic
    #fig = plot(aspect_ratio = 1,legend = false)
    #AS.plot_heuristic!(heuristic)
    #display(fig)
end

function BB.compute_lower_bound!(prob::OptimalControlProblem, node::BB.Node)
    path = node.elem
    cells = prob.cells
    curr_cell = cells[path[end]]
    if length(path) == 1
        node.lower_bound = curr_cell.lower_bound
        return
    end
    prev_cell = cells[path[end-1]]
    from = length(path) == 2 ? -1 : path[end-2]
    #if the alternating simulation is not yet computed
    if prev_cell.coarse_abstraction==nothing
        build_coarse_abstraction!(prob,prev_cell)
    end
    # build the heuristic if necessary
    if !haskey(prev_cell.heuristics,from)
        build_heuristic!(prev_cell,from)
    end
    if curr_cell.index == prob.qT
        if curr_cell.coarse_abstraction==nothing
            build_coarse_abstraction!(prob,curr_cell)
        end
        if !haskey(curr_cell.heuristics,path[end-1])
            build_heuristic!(curr_cell,path[end-1])
        end
    end
    # compute the lower bound
    heuristic = prev_cell.heuristics[from]
    _T_ = prev_cell.local_target_set[path[end]]
    a = AS.get_min_value_heurisitic(heuristic,_T_)
    cost = node.parent.lower_bound + a
    if prev_cell.index != prob.qT
        cost -= prev_cell.lower_bound
    else
        heuristic = prev_cell.heuristics[path[end-2]]
        _T_ = prev_cell.local_target_set[-1]
        b = AS.get_min_value_heurisitic(heuristic,_T_)
        cost -= b
    end
    if curr_cell.index != prob.qT
        cost += curr_cell.lower_bound
    else
        heuristic = curr_cell.heuristics[path[end-1]]
        _T_ = curr_cell.local_target_set[-1]
        b = AS.get_min_value_heurisitic(heuristic,_T_)
        cost += b
    end
    #println(a)
    #println("cost:  ",path," -> ",cost)
    node.lower_bound = cost + 2.0#+ 5 # where 1 is the cost of the transition
end

function BB.compute_upper_bound!(prob::OptimalControlProblem, node::BB.Node)
    path = node.elem
    cells = prob.cells
    if path[end]!=prob.qT
        return
    end
    for (i,index) in enumerate(path)
        cell = cells[index]
        # if the abstraction (or any part of it) is not yet build
        if cell.fine_abstraction == nothing
            hx = SVector(0.5, 0.5)*1.0; x0 = @SVector zeros(prob.dim)
            Xgrid = AB.GridFree(x0,hx)
            Xdom = AB.DomainList(Xgrid)
            AB.add_set!(Xdom, cells[index].hyperrectangle , AB.INNER)
            for i in cell.outneighbors
                AB.add_set!(Xdom,cell.local_target_set[i], AB.INNER)
            end
            symmodel = AB.NewSymbolicModelListList(Xdom, prob.Udom)
            cell.fine_abstraction = (symmodel,nothing)
        end

        prev = i==1 ? -1 : path[i-1]
        next = i==length(path) ? -1 : path[i+1]
        controllers = cells[index].controllers
        # build the controller if not already computed
        if !haskey(controllers,(prev,next))
            (symmodel, transitions_previously_added) = cell.fine_abstraction
            _I_ = cell.local_init_set[prev]
            initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
            _T_ = cell.local_target_set[next]
            targetlist = U.get_symbol(symmodel,_T_,AB.INNER)
            heuristic = cell.heuristics[prev]

            problem = LA.compute_controller(symmodel, prob.contsys, initlist, targetlist,prob.transition_cost, prob.functions[4], prob.functions[3], h2, heuristic_data=heuristic,transitions_previously_added=transitions_previously_added)
            cell.fine_abstraction = (symmodel,problem.transitions_previously_added)
            contr = problem.contr
            controllers[(prev,next)] = contr
            #LA.plot_result!(problem)
        end
    end
    # compute the upper bound
    (traj, cost, succeed) = simulate_trajectory(prob, path)
    node.upper_bound = cost
    node.sol = path
    println()
    println(cost)
    fig = plot(aspect_ratio = 1,legend = false)
    for cell in prob.cells
        h = cell.hyperrectangle
        plot!(U.rectangle(h.lb,h.ub), opacity=.4,color=:blue)
    end
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
        if prob.option && (f1(reverse(new_path)) || (length(path) >= 2 && neighbor_cell == path[end-1]) || f2(new_path))
            continue
        end
        push!(children,BB.Node(new_path,node.depth+1,parent=node,lower_bound=-Inf, upper_bound=Inf,ext=prob)) #prob : temporary (to change)
    end
    return children
end

function simulate_trajectory(prob::OptimalControlProblem, path)
    x0 = prob.x0
    traj = []
    cost = 0.0

    for (i,index) in enumerate(path)
        cell = prob.cells[index]
        (symmodel,tab) = cell.fine_abstraction
        prev = i==1 ? -1 : path[i-1]
        next = i==length(path) ? -1 : path[i+1]
        contr = cell.controllers[(prev,next)]

        _T_ = cell.local_target_set[next]
        targetlist = U.get_symbol(symmodel,_T_,AB.INNER)
        (local_traj,local_cost,succeed) = simulate_local_trajectory(prob.contsys, symmodel, contr, x0, targetlist; transition_cost=prob.transition_cost)
        if !succeed
            return (traj,cost,false)
        end
        traj = [traj..., local_traj[2:end]...]
        cost += local_cost
        x0 = local_traj[end]
    end
    return (traj,cost,true)
end
function simulate_local_trajectory(contsys, symmodel, contr, x0, targetlist; transition_cost)
    traj = []
    cost = 0.0
    while true
        push!(traj,x0)
        xpos = AB.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom) # note, this should normally never happen
            @warn("Trajectory out of domain")
            return (traj,cost,false)
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        if source ∈ targetlist
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
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return (traj,cost,true)
end

function print_trajectory!(traj)
    for i=1:length(traj)-1
        plot!([traj[i][1],traj[i+1][1]], [traj[i][2],traj[i+1][2]],color =:red,linewidth = 2)
        if i>1
            scatter!([traj[i][1]],[traj[i][2]],color =:red,markersize=2)
        end
    end
    scatter!([traj[1][1]],[traj[1][2]],color =:green,markersize=3)
    scatter!([traj[end][1]],[traj[end][2]],color =:yellow,markersize=3)
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
