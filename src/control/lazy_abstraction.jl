module LazyAbstractionReach
"""
Module to solve a reachability problem.
For this, it builds an abstraction that is alternatingly simulated by our original system.
This abstraction and the controller are built lazily based on a heuristic.
"""


using ..Utils
UT = Utils

using ..Domain
DO = Domain

using ..Symbolic
SY = Symbolic

using ..Control
CO = Control

using StaticArrays, Plots

struct State
    source::Int
end
Base.:(==)(s1::State, s2::State)= s1.source == s2.source

mutable struct MutableMatrix{T, VT<:AbstractVector{T}}
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
function Base.setindex!(m::MutableMatrix{T}, val::T, col::Int, row::Int) where T
    m.data[(col - 1) * m.num_rows + row] = val
end

mutable struct LazyAbstraction{T,SM,C,TC<:Function,PrI<:Function,PoI<:Function,HD} <: UT.SearchProblem{T}
    initial::Vector{T} # list of initial states of A* (target set)
    goal::Vector{T}    # list of goal states of A* (init set)
    symmodel::SM
    contsys::C
    transition_cost::TC #transition_cost(x,u) function
    pre_image::PrI  # function to compute the list of potential pre-image of a cell for a given input
    post_image::PoI # function to compute the list of potential post-image of a cell for a given input
    transitions_added::MutableMatrix{Bool,BitVector}          # could be an array or a dictionnary (to be added)
    num_targets_unreachable::MutableMatrix{Int,Vector{Int}}     # could be an array or a dictionnary (to be added)
    controllable::BitVector                # could be an array or a dictionnary (to be added)
    num_init_unreachable::Int   # counter of the remaining non controllable init cells
    heuristic_data::HD # extension for potential additionnal data for the heuristic function
    contr::UT.SortedTupleSet{2,Int}  # controller
    closed::Union{Nothing,Dict{T,Bool}} # only usefull for the printing (could be discard later)
    costs_temp::MutableMatrix{Float64,Vector{Float64}} # array containing the current worse cost to reach the target, if the next input applied is symbol
    costs::Vector{Float64} # vector containing the (worst) cost to reach the target set for each cell (necessary because of the pseudo non determinism) = Lyapunov function
    transitions_previously_added::MutableMatrix{Int,Vector{Int}} # only necessary, if we need to reuse a partially computed symmodel
                                                                 # this array should contain the number of outgoing neighbors for
                                                                 # previously computed couple (cell,input) and -1 for the others.
end

function LazyAbstraction(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image;heuristic_data=nothing,transitions_previously_added=nothing)
    transitions_added = MutableMatrix(false, symmodel.autom.nsymbols, symmodel.autom.nstates)
    num_targets_unreachable = MutableMatrix(0, symmodel.autom.nsymbols, symmodel.autom.nstates)
    controllable = falses(symmodel.autom.nstates)
    for init in initial
        controllable[init.source] = true
    end
    num_init_unreachable = length(goal)
    costs_temp = MutableMatrix(0.0, symmodel.autom.nsymbols, symmodel.autom.nstates)
    costs = zeros(Float64, symmodel.autom.nstates)
    contr =  CO.NewControllerList()
    if transitions_previously_added === nothing
        transitions_previously_added = MutableMatrix(-1, symmodel.autom.nsymbols, symmodel.autom.nstates)
    end
    return LazyAbstraction{State,typeof(symmodel),typeof(contsys),typeof(transition_cost),typeof(pre_image),typeof(post_image),typeof(heuristic_data)}(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image,transitions_added,num_targets_unreachable,controllable,num_init_unreachable,heuristic_data,contr,nothing,costs_temp,costs,transitions_previously_added)
end

function UT.goal_test(problem::LazyAbstraction, state::State)
    if state.source in [s.source for s in problem.goal]
        if iszero(problem.num_init_unreachable-=1)
            return true
        end
    end
    return false
end

# for the moment, one action costs 1
function UT.path_cost(problem::LazyAbstraction,c, state1::State, action, state2::State)
    source = state2.source
    pos = SY.get_xpos_by_state(problem.symmodel,source)
    x = DO.get_coord_by_pos(problem.symmodel.Xdom.grid,pos)
    upos = SY.get_upos_by_symbol(problem.symmodel,action)
    u = DO.get_coord_by_pos(problem.symmodel.Udom.grid,upos)

    problem.costs[source] += problem.transition_cost(x,u) #add the cost of the transition (should be the worst for the cell)
    return problem.costs[source] #c + 0.5
end

function transitions!(source,symbol,u,symmodel,contsys,post_image)
    xpos = SY.get_xpos_by_state(symmodel, source)
    over_approx = post_image(symmodel,contsys,xpos,u)
    translist = [(cell, source, symbol) for cell in over_approx]
    SY.add_transitions!(symmodel.autom, translist)
    return length(over_approx)
end

function _update_cache!(problem, ns1, ns2, nsym)
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

function update_abstraction!(successors,problem,source)
    symmodel = problem.symmodel
    contsys = problem.contsys
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    nsym = symmodel.autom.nsymbols

    xpos = SY.get_xpos_by_state(symmodel, source)
    x = DO.get_coord_by_pos(Xdom.grid, xpos)
    for upos in DO.enum_pos(Udom)
        symbol = SY.get_symbol_by_upos(symmodel, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)

        ns1 = symmodel.autom.nstates
        cells = problem.pre_image(symmodel,contsys,xpos,u)
        _update_cache!(problem, ns1, symmodel.autom.nstates, nsym)
        for cell in cells
            if !problem.controllable[cell]
                if !problem.transitions_added[cell,symbol]
                    # add transitions for input u starting from cell if it has not be done yet
                    n_trans = 0
                    if problem.transitions_previously_added[cell,symbol] != -1
                        n_trans = problem.transitions_previously_added[cell,symbol]
                    else
                        ns1 = symmodel.autom.nstates
                        n_trans = transitions!(cell,symbol,u,symmodel,contsys,problem.post_image)
                        _update_cache!(problem, ns1, symmodel.autom.nstates, nsym)
                        problem.transitions_previously_added[cell,symbol] = n_trans
                    end
                    problem.num_targets_unreachable[cell,symbol] = n_trans
                    problem.transitions_added[cell,symbol] = true
                end
                # check if the cell is really in the pre-image
                if (source,cell,symbol) in symmodel.autom.transitions
                    problem.costs_temp[cell,symbol] = max(problem.costs_temp[cell,symbol],problem.costs[source])
                    if iszero(problem.num_targets_unreachable[cell,symbol] -= 1)
                        problem.costs[cell] = problem.costs_temp[cell,symbol]
                        problem.controllable[cell] = true
                        push!(successors,(symbol,State(cell)))
                        UT.push_new!(problem.contr, (cell, symbol))
                    end
                end
            end
        end
    end
end

function UT.successor(problem::LazyAbstraction, state::State)
    successors = []
    update_abstraction!(successors,problem,state.source)
    return successors
end

function compute_controller(symmodel, contsys, initlist::Vector{Int}, targetlist::Vector{Int}, transition_cost, pre_image, post_image, h; heuristic_data=nothing,transitions_previously_added=nothing)
    initial = [State(tar) for tar in targetlist]
    goal = [State(init) for init in initlist]
    problem = LazyAbstraction(initial,goal,symmodel,contsys,transition_cost,pre_image,post_image,heuristic_data=heuristic_data,transitions_previously_added=transitions_previously_added)
    node, nb = UT.astar_graph_search(problem,h)
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
        if any(u -> problem.transitions_added[k,u], 1:problem.transitions_added.num_rows)
            pos = SY.get_xpos_by_state(symmodel, k)
            if !haskey(dict,pos[dims])
                dict[pos[dims]] = true
                center = DO.get_coord_by_pos(grid, pos)
                Plots.plot!(rectangle(center[dims],h./2), opacity=.2,color=:yellow)
            end
        end
    end

    # controllable state
    dict = Dict{NTuple{2,Int}, Any}()
    for (cell, symbol) in contr.data
        pos = SY.get_xpos_by_state(symmodel,cell)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims],h./2), opacity=.3,color=:blue)
        end
    end

    # states selected by A* to compute their pre-image
    dict = Dict{NTuple{2,Int}, Any}()
    for state in Base.keys(problem.closed)
        pos = SY.get_xpos_by_state(symmodel,state.source)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims],h./2), opacity=.5,color=:blue)
        end
    end

    # initial set
    dict = Dict{NTuple{2,Int}, Any}()
    for s in initlist
        pos = SY.get_xpos_by_state(symmodel,s)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims],h./2), opacity=.4,color=:green)
        end
    end

    # target set
    dict = Dict{NTuple{2,Int}, Any}()
    for s in targetlist
        pos = SY.get_xpos_by_state(symmodel,s)
        if !haskey(dict,pos[dims])
            dict[pos[dims]] = true
            center = DO.get_coord_by_pos(grid, pos)
            Plots.plot!(rectangle(center[dims],h./2), opacity=.5,color=:red)
        end
    end

    # plot a trajectory
    if x0 != nothing
        (traj,success) = trajectory_reach(contsys, symmodel, contr, x0, targetlist)
        plot_trajectory!(symmodel,traj,dims=dims)
    end
end


function trajectory_reach(contsys, symmodel, contr, x0, targetlist; randchoose = false)
    traj = []
    while true
        x0 = DO.set_in_period_coord(symmodel.Xdom,x0)
        push!(traj,x0)
        xpos = DO.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return (traj,false)
        end
        source = SY.get_state_by_xpos(symmodel, xpos)[1]
        if source ∈ targetlist
            break
        end
        symbollist = UT.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return (traj,false)
        end
        if randchoose
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end

        upos = SY.get_upos_by_symbol(symmodel, symbol)
        u = DO.get_coord_by_pos(symmodel.Udom.grid, upos)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return (traj,true)
end


function plot_trajectory!(symmodel,traj;dims=[1,2])
    domain = symmodel.Xdom
    grid = domain.grid
    k = dims[1]; l = dims[2]
    for i=1:length(traj)-1
        Plots.plot!([traj[i][k],traj[i+1][k]], [traj[i][l],traj[i+1][l]],color =:red,linewidth = 2)
        if i>1
            Plots.scatter!([traj[i][k]],[traj[i][l]],color =:red,markersize=2)
        end
    end
    Plots.scatter!([traj[1][k]],[traj[1][l]],color =:green,markersize=3)
    Plots.scatter!([traj[end][k]],[traj[end][l]],color =:yellow,markersize=3)
end

end #end module
