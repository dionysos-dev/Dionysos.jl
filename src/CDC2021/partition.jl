module Partition

using ..Abstraction
AB = Abstraction

using ..DomainList
D = DomainList

using ..Utils
U = Utils

using ..AlternatingSimulation
AS = AlternatingSimulation

using StaticArrays, LightGraphs, Plots

mutable struct Cell
    index::Int
    hyperrectangle::AB.HyperRectangle
    outneighbors # list of out neighbors
    inneighbors

    medium_abstraction # symmmodel of the alternating system
    fine_abstraction   # symmmodel of the alternatingly simulated system

    heuristics::Dict{Int,Any} # only depends of the previous cell
    controllers::Dict{NTuple{2,Int},Any} # depends from the previous and the next cells

    lower_bound::Float64 # cost of the big cell abstraction
    upper_bound::Float64 # (probably not necessary, not used)

    reachable_set
    local_target_set
    local_init_set
end

function Initialise(partition,contsys,Udom,_I_,_T_,compute_reachable_set,minimum_transition_cost,periodic,periods)
    X,hx = partition
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, X , AB.INNER)
    U.plot_domain!(Xdom,dims=[1,2],opacity=0.15,color=:blue)
    symmodel = AB.NewSymbolicModelListList(Xdom,Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_2)
    symmodel.autom = AS.build_alternating_simulation(problem)

    q0 = localise(symmodel,_I_)
    qT = localise(symmodel,_T_)

    heuristic = AS.build_heuristic(symmodel,[qT]) #remark later, we could put the several which contains the targetset
    #AS.plot_heuristic!(heuristic,opacity=0.15,color=:blue)

    grid = Xdom.grid
    r = grid.h/2
    cells = Cell[]
    for pos in AB.enum_pos(Xdom)
        i = AB.get_state_by_xpos(symmodel,pos)
        outneighbors = LightGraphs.outneighbors(symmodel.autom,i)
        inneighbors = LightGraphs.inneighbors(symmodel.autom,i)
        n = length(outneighbors)
        medium_abstraction =  nothing
        fine_abstraction = nothing
        heuristics = Dict{Int, Any}()
        controllers = Dict{NTuple{2,Int}, Any}()
        lower_bound = heuristic.dists[i]
        upper_bound = Inf
        c = AB.get_coord_by_pos(grid, pos)
        rec = AB.HyperRectangle(c-r,c+r)
        reachable_set = compute_reachable_set(rec,contsys,Udom)
        local_target_set = Dict{Int, Any}()
        local_init_set = Dict{Int, Any}()
        cell = Cell(i,rec,outneighbors,inneighbors,medium_abstraction,fine_abstraction,heuristics,controllers,
                    lower_bound,upper_bound,reachable_set,local_target_set,local_init_set)
        push!(cells,cell)
    end
    compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set,periodic,periods)

    plot_local_set(cells[1],Xdom)#587

    return (q0,qT,symmodel,cells)
end

function localise(symmodel,_I_)
    Xdom = symmodel.Xdom
    rectI = AB.get_pos_lims(Xdom.grid, _I_, AB.OUTER)
    pos_iter = Iterators.product(AB._ranges(rectI)...)
    q = []
    for pos in pos_iter
        pos = D.set_in_period_pos(Xdom,pos)
        if pos in Xdom
            push!(q,AB.get_state_by_xpos(symmodel,pos))
        end
    end
    return q[1]
end

function one_direction(lb,ub,T)
    if ub-lb>=T
        return [[0,T]]
    else
        lb = mod(lb,T)
        ub = mod(ub,T)
        if lb<=ub
            return [[lb,ub]]
        else
            return [[0,ub],[lb,T]]
        end
    end
end
function recursive(L,rec,lb,ub,periodic,periods,i)
    if i>length(periodic)
        N = length(lb)
        push!(L,AB.HyperRectangle(SVector{N}(lb), SVector{N}(ub)))
        return
    end
    dim = periodic[i]
    intervals = one_direction(rec.lb[dim],rec.ub[dim],periods[i])
    for interval in intervals
        l = copy(lb)
        u = copy(ub)
        l[dim] = interval[1]
        u[dim] = interval[2]
        recursive(L,rec,l,u,periodic,periods,i+1)
    end
end
function set_rec_in_period(periodic,periods,rec)
    L = []
    lb = collect(rec.lb)
    ub = collect(rec.ub)
    recursive(L,rec,lb,ub,periodic,periods,1)
    return L
end

function compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set,periodic,periods)
    cells[q0].local_init_set[-1] = []
    cells[qT].local_target_set[-1] = []
    for cell in cells
        for i in cell.inneighbors
            cell.local_init_set[i] = []
        end
        for i in cell.outneighbors
            cell.local_target_set[i] = []
        end
    end
    push!(cells[q0].local_init_set[-1],set_rec_in_period(periodic,periods,_I_)...)
    push!(cells[qT].local_target_set[-1],set_rec_in_period(periodic,periods,_T_)...)
    for cell in cells
        R = compute_reachable_set(cell.hyperrectangle,contsys,Udom)
        RL = set_rec_in_period(periodic,periods,R)
        for rec in RL
            for i in cell.outneighbors
                neigh = cells[i]
                I = AB.intersect(rec, neigh.hyperrectangle)
                push!(neigh.local_init_set[cell.index],I)
                push!(cell.local_target_set[neigh.index],I)
            end
        end
    end
end

function plot_local_set(cell,Xdom;dims=[1,2])
    rec = cell.hyperrectangle
    reachable_set = cell.reachable_set
    reachable_set = set_rec_in_period(Xdom.periodic,Xdom.periods,reachable_set)

    fig = plot(aspect_ratio = 1,legend = false)
    dims=[1,2]
    U.plot_domain!(Xdom,dims=dims)
    plot!(U.rectangle(rec.lb[dims],rec.ub[dims]), opacity=0.2,color=:red)
    for rect in reachable_set
        plot!(U.rectangle(rect.lb[dims],rect.ub[dims]), opacity=0.4,color=:red)
    end
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    dims=[3,4]
    U.plot_domain!(Xdom,dims=dims)
    plot!(U.rectangle(rec.lb[dims],rec.ub[dims]), opacity=0.2,color=:red)
    for rect in reachable_set
        plot!(U.rectangle(rect.lb[dims],rect.ub[dims]), opacity=0.4,color=:red)
    end
    display(fig)
end
end # end module
