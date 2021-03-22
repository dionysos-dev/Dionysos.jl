module Lazy_abstraction
"""
Module that aims to solve a reachability problem.
For this, he builds an abstraction that is alternatingly simulated by our original system.
This abstraction and the controller are built lazily based on a heuristic.
(Similar to the module Abstraction except that the computation are done lazily).
"""
#TO DO: -array vs dictionnary
#       -eventually delete State struct

include("search.jl")
using .Search
S = Search

using ..Abstraction
AB = Abstraction

using ..DomainList
D = DomainList

using Plots,StaticArrays

struct State
    source
end
Base.:(==)(s1::State, s2::State)= s1.source == s2.source

mutable struct LazyAbstraction{T} <: S.SearchProblem{T}
    initial # list of initial states of A* (target set)
    goal    # list of goal states of A* (init set)
    symmodel
    contsys
    transition_cost #transition_cost(x,u) function
    pre_image  # function to compute the list of potential pre-image of a cell for a given input
    post_image # function to compute the list of potential post-image of a cell for a given input
    transitions_added::Union{Array{Bool,2}, DataType}          # could be an array or a dictionnary (to be added)
    num_targets_unreachable::Union{Array{Int,2}, DataType}     # could be an array or a dictionnary (to be added)
    controllable::Union{Vector{Bool},DataType}                 # could be an array or a dictionnary (to be added)
    num_init_unreachable::Int   # counter of the remaining non controllable init cells
    heuristic_data # extension for potential additionnal data for the heuristic function
    contr  # controller
    closed # only usefull for the printing (could be discard later)
    costs_temp::Union{Array{Float64,2}, DataType} # array containing the current worse cost to reach the target, if the next input applied is symbol
    costs::Union{Vector{Float64}, DataType} # vector containing the (worst) cost to reach the target set for each cell (necessary because of the pseudo non determinism) = Lyapunov function
    transitions_previously_added::Union{Array{Int,2}, DataType} # only necessary, if we need to reuse a partially computed symmodel
                                                                # this array should contain the number of outgoing neighbors for
                                                                # previously computed couple (cell,input) and -1 for the others.
end

function LazyAbstraction(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image;heuristic_data=nothing,transitions_previously_added=nothing)
    transitions_added = fill(false, (symmodel.autom.nstates,symmodel.autom.nsymbols))
    num_targets_unreachable = zeros(Int, symmodel.autom.nstates, symmodel.autom.nsymbols)
    controllable = fill(false, (symmodel.autom.nstates))
    for init in initial
        controllable[init.source] = true
    end
    num_init_unreachable = length(goal)
    costs_temp = zeros(Float64, symmodel.autom.nstates, symmodel.autom.nsymbols)
    costs = zeros(Float64, symmodel.autom.nstates)
    contr =  AB.NewControllerList()
    if transitions_previously_added == nothing
        transitions_previously_added = fill(-1, (symmodel.autom.nstates,symmodel.autom.nsymbols))
    end
    return LazyAbstraction{State}(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image,transitions_added,num_targets_unreachable,controllable,num_init_unreachable,heuristic_data,contr,nothing,costs_temp,costs,transitions_previously_added)
end

function S.goal_test(problem::LazyAbstraction, state::State)
    if state.source in [s.source for s in problem.goal]
        if iszero(problem.num_init_unreachable-=1)
            return true
        end
    end
    return false
end

# for the moment, one action costs 1
function S.path_cost(problem::LazyAbstraction,c, state1::State, action, state2::State)
    source = state2.source
    pos = AB.get_xpos_by_state(problem.symmodel,source)
    x = AB.get_coord_by_pos(problem.symmodel.Xdom.grid,pos)
    upos = AB.get_upos_by_symbol(problem.symmodel,action)
    u = AB.get_coord_by_pos(problem.symmodel.Udom.grid,upos)

    problem.costs[source] += problem.transition_cost(x,u) #add the cost of the transition (should be the worst for the cell)
    return problem.costs[source] #c + 0.5
end

function transitions!(source,symbol,u,symmodel,contsys,post_image)
    xpos = AB.get_xpos_by_state(symmodel, source)
    over_approx = post_image(symmodel,contsys,xpos,u)
    translist = [(cell, source, symbol)  for cell in over_approx]
    AB.add_transitions!(symmodel.autom, translist)
    return length(over_approx)
end

function update_abstraction!(successors,problem,source)
    symmodel = problem.symmodel
    contsys = problem.contsys
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom

    xpos = AB.get_xpos_by_state(symmodel, source)
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    for upos in AB.enum_pos(Udom)
        symbol = AB.get_symbol_by_upos(symmodel, upos)
        u = AB.get_coord_by_pos(Udom.grid, upos)

        cells = problem.pre_image(symmodel,contsys,xpos,u)
        for cell in cells
            if !problem.controllable[cell]
                if !problem.transitions_added[cell,symbol]
                    # add transitions for input u starting from cell if it has not be done yet
                    n_trans = 0
                    if problem.transitions_previously_added[cell,symbol] != -1
                        n_trans = problem.transitions_previously_added[cell,symbol]
                    else
                        n_trans = transitions!(cell,symbol,u,symmodel,contsys,problem.post_image)
                        problem.transitions_previously_added[cell,symbol] = n_trans
                    end
                    problem.num_targets_unreachable[cell,symbol] = n_trans
                    problem.transitions_added[cell,symbol] = true
                end
                # check if the cell is really in the pre-image
                if (source,cell,symbol) in symmodel.autom.transitions.data
                    #println("in the pre-image")
                    problem.costs_temp[cell,symbol] = max(problem.costs_temp[cell,symbol],problem.costs[source])
                    if iszero(problem.num_targets_unreachable[cell,symbol] -= 1)
                        println("cell added (controlled)")
                        problem.costs[cell] = problem.costs_temp[cell,symbol]
                        problem.controllable[cell] = true
                        push!(successors,(symbol,State(cell)))
                        AB.push_new!(problem.contr, (cell, symbol))
                    end
                end
            end
        end
    end
end

function S.successor(problem::LazyAbstraction, state::State)
    successors = []
    #readline()
    #fig = plot(aspect_ratio = 1,legend = false)
    #plot_result!(problem,dims=[1,2])
    #display(fig)
    update_abstraction!(successors,problem,state.source)
    return successors
end

function compute_controller(symmodel, contsys, initlist::Vector{Int}, targetlist::Vector{Int}, transition_cost, pre_image, post_image, h; heuristic_data=nothing,transitions_previously_added=nothing)
    initial = [State(tar) for tar in targetlist]
    goal = [State(init) for init in initlist]
    problem = LazyAbstraction(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image,heuristic_data=heuristic_data,transitions_previously_added=transitions_previously_added)
    node, nb = S.astar_graph_search(problem,h)
    println("\nnumber of transitions created: ", length(problem.symmodel.autom.transitions))
    if node == nothing
        println("compute_controller_reach! terminated without covering init set")
        return problem,false
    end
    println("compute_controller_reach! terminated with success")
    return problem,true
end



## printing
function rectangle(c,r)
    Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
end

function plot_result!(problem;dims=[1,2],x0=nothing)
    println()
    println("Plotting")
    targetlist = [init.source for init in problem.initial]
    initlist = [goal.source for goal in problem.goal]
    contsys = problem.contsys
    contr = problem.contr
    symmodel = problem.symmodel
    domain = symmodel.Xdom
    grid = domain.grid
    h = grid.h[dims]

    # states for which transisitons have been computed for at least one input
    dict = Dict{NTuple{2,Int}, Any}()
    for k = 1:symmodel.autom.nstates
        if any(problem.transitions_added[k,:])
            pos = AB.get_xpos_by_state(symmodel, k)
            if !haskey(dict,pos[dims])
                dict[pos[dims]] = true
                center = AB.get_coord_by_pos(grid, pos)
                plot!(rectangle(center[dims],h./2), opacity=.2,color=:yellow)
            end
        end
    end

    # controllable state
    dict = Dict{NTuple{2,Int}, Any}()
    for (cell, symbol) in contr.data
        pos = AB.get_xpos_by_state(symmodel,cell)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = AB.get_coord_by_pos(grid, pos)
            plot!(rectangle(center[dims],h./2), opacity=.3,color=:blue)
        end
    end

    # states selected by A* to compute their pre-image
    dict = Dict{NTuple{2,Int}, Any}()
    for state in Base.keys(problem.closed)
        pos = AB.get_xpos_by_state(symmodel,state.source)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = AB.get_coord_by_pos(grid, pos)
            plot!(rectangle(center[dims],h./2), opacity=.5,color=:blue)
        end
    end

    # initial set
    dict = Dict{NTuple{2,Int}, Any}()
    for s in initlist
        pos = AB.get_xpos_by_state(symmodel,s)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = AB.get_coord_by_pos(grid, pos)
            plot!(rectangle(center[dims],h./2), opacity=.4,color=:green)
        end
    end

    # target set
    dict = Dict{NTuple{2,Int}, Any}()
    for s in targetlist
        pos = AB.get_xpos_by_state(symmodel,s)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = AB.get_coord_by_pos(grid, pos)
            plot!(rectangle(center[dims],h./2), opacity=.5,color=:red)
        end
    end

    # plot a trajectory
    if x0 != nothing
        (traj,success) = trajectory_reach(contsys, symmodel, contr, x0, targetlist)
        print_trajectory!(symmodel,traj,dims=dims)
    end
end


function trajectory_reach(contsys, symmodel, contr, x0, targetlist; randchoose = false)
    traj = []
    while true
        x0 = D.set_in_period_coord(symmodel.Xdom,x0)
        push!(traj,x0)
        xpos = AB.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return (traj,false)
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        if source ∈ targetlist
            break
        end
        symbollist = AB.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return (traj,false)
        end
        if randchoose
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end

        upos = AB.get_upos_by_symbol(symmodel, symbol)
        u = AB.get_coord_by_pos(symmodel.Udom.grid, upos)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return (traj,true)
end


function print_trajectory!(symmodel,traj;dims=[1,2])
    domain = symmodel.Xdom
    grid = domain.grid
    k = dims[1]; l = dims[2]
    for i=1:length(traj)-1
        plot!([traj[i][k],traj[i+1][k]], [traj[i][l],traj[i+1][l]],color =:red,linewidth = 2)
        if i>1
            scatter!([traj[i][k]],[traj[i][l]],color =:red,markersize=2)
        end
    end
    scatter!([traj[1][k]],[traj[1][l]],color =:green,markersize=3)
    scatter!([traj[end][k]],[traj[end][l]],color =:yellow,markersize=3)
end
end # end module
