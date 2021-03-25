module Partition

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const U = Utils

using ..AlternatingSimulation
const AS = AlternatingSimulation

using StaticArrays, LightGraphs, Plots

mutable struct Cell{N,T}
    index::Int
    hyperrectangle::AB.HyperRectangle
    outneighbors::SubArray{Int,1,Array{Int,1},Tuple{UnitRange{Int}},true} # list of out neighbors
    inneighbors::Vector{Int}

    medium_abstraction # symmmodel of the alternating system
    fine_abstraction   # symmmodel of the alternatingly simulated system

    heuristics::Dict{Int,Any} # only depends of the previous cell
    controllers::Dict{NTuple{2,Int},Any} # depends from the previous and the next cells

    lower_bound::Float64 # cost of the big cell abstraction
    upper_bound::Float64 # (probably not necessary, not used)

    reachable_set::AB.HyperRectangle{SVector{N,T}}
    local_target_set::Dict{Int,Vector{AB.HyperRectangle{SVector{N,T}}}}
    local_init_set::Dict{Int,Vector{AB.HyperRectangle{SVector{N,T}}}}
end


function compute_reachable_sets(Xdom::AB.Domain{N,T},contsys,Udom,compute_reachable_set) where {N,T}
    grid = Xdom.grid
    L = AB.HyperRectangle{SVector{N,T}}[]
    i = 1
    nc = AB.get_ncells(Xdom)
    h = grid.h
    for pos in AB.enum_pos(Xdom)
        println(i," / ",nc)
        c = AB.get_coord_by_pos(grid,pos)
        hyperrectangle = AB.HyperRectangle(c-h/2,c+h/2)
        reachable_set = compute_reachable_set(hyperrectangle,contsys,Udom)
        push!(L,reachable_set)
        i+=1
    end
    return L
end

function Initialise(partition,contsys,Udom,_I_,_T_,compute_reachable_set,minimum_transition_cost,periodic,periods,T0)
    X,hx,O = partition

    grid = D.build_grid_in_rec(X,hx)
    Xdom = D.GeneralDomainList(grid;periodic=periodic,periods=periods,T0=T0,lims=X)
    #Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, X , AB.INNER)
    for obstacle in O
        AB.remove_set!(Xdom, obstacle, AB.OUTER)
    end
    fig = plot(aspect_ratio = 1,legend = false)
    U.plot_domain!(Xdom,dims=[1,2],opacity=0.15,color=:blue)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    U.plot_domain!(Xdom,dims=[3,4],opacity=0.15,color=:blue)
    display(fig)
    L = compute_reachable_sets(Xdom,contsys,Udom,compute_reachable_set)

    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_3,ext=L)
    autom = AS.build_alternating_simulation(problem)
    symmodel = AB.with_automaton(symmodel, autom)

    q0 = localise(symmodel,_I_)
    qT = localise(symmodel,_T_)
    println(q0)
    println(qT)
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
        medium_abstraction = nothing
        fine_abstraction = nothing
        heuristics = Dict{Int, Any}()
        controllers = Dict{NTuple{2,Int}, Any}()
        lower_bound = heuristic.dists[i]
        upper_bound = Inf
        c = AB.get_coord_by_pos(grid, pos)
        rec = AB.HyperRectangle(c-r,c+r)
        reachable_set = L[i]
        local_target_set = Dict{Int, Vector{typeof(_I_)}}()
        local_init_set = Dict{Int, Vector{typeof(_I_)}}()
        cell = Cell(i,rec,outneighbors,inneighbors,medium_abstraction,fine_abstraction,heuristics,controllers,
                    lower_bound,upper_bound,reachable_set,local_target_set,local_init_set)
        push!(cells,cell)
    end
    compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set,periodic,periods,T0)

    for i = 1:15
        plot_local_set(cells[i],Xdom)#587
        plot_local_set(cells[i],Xdom;dims=[3,4])
    end
    return (q0,qT,symmodel,cells)
    #return 1,1,1,1
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

function compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set,periodic,periods,T0)
    cells[q0].local_init_set[-2] = []
    cells[qT].local_target_set[-1] = []
    for cell in cells
        for i in cell.inneighbors
            cell.local_init_set[i] = []
        end
        for i in cell.outneighbors
            cell.local_target_set[i] = []
        end
    end
    push!(cells[q0].local_init_set[-2],D.set_rec_in_period(periodic,periods,T0,_I_)...)
    push!(cells[qT].local_target_set[-1],D.set_rec_in_period(periodic,periods,T0,_T_)...)
    for cell in cells
        RL = D.set_rec_in_period(periodic,periods,T0,cell.reachable_set)
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
    reachable_set = D.set_rec_in_period(Xdom.periodic,Xdom.periods,Xdom.T0,reachable_set)

    fig = plot(aspect_ratio = 1,legend = false)
    U.plot_domain!(Xdom,dims=dims)
    plot!(U.rectangle(rec.lb[dims],rec.ub[dims]), opacity=0.2,color=:red)
    for rect in reachable_set
        plot!(U.rectangle(rect.lb[dims],rect.ub[dims]), opacity=0.4,color=:red)
    end
    display(fig)
end
end # end module
