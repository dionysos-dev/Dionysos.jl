module Partition

using ..Abstraction
AB = Abstraction

using ..Utils
U = Utils

using ..AlternatingSimulation
AS = AlternatingSimulation

using StaticArrays, LightGraphs, Plots

mutable struct Cell
    index::Int
    hyperrectangle::AB.HyperRectangle
    outneighbors # list of out neighbors

    coarse_abstraction # symmmodel of the alternating system
    fine_abstraction   # symmmodel of the alternatingly simulated system

    heuristics::Dict{Int,Any} # only depends of the previous cell
    controllers::Dict{NTuple{2,Int},Any} # depends from the previous and the next cells

    lower_bound::Float64 # cost of the big cell abstraction
    upper_bound::Float64 # (probably not necessary, not used)

    local_target_set
    local_init_set
end

function localise(rectangles::Vector{AB.HyperRectangle{T}},_T_::AB.HyperRectangle{T}) where T
    for (i,rec) in enumerate(rectangles)
        if U.is_intersection(rec, _T_)
            return i
        end
    end
end
function Initialise(rectangles,contsys,Udom,_I_,_T_,compute_reachable_set,minimum_transition_cost)
    q0 = localise(rectangles,_I_)
    qT = localise(rectangles,_T_)
    cells = build_coarser_abstraction(rectangles,contsys,Udom,qT,compute_reachable_set,minimum_transition_cost)
    compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set)
    return (q0,qT,cells)
end

function build_coarser_abstraction(rectangles,contsys,Udom,qT,compute_reachable_set,minimum_transition_cost)
    Xdom = U.CustomList(rectangles)
    symmodel = U.SymbolicModelList2(Xdom,Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_1)
    symmodel.autom = AS.build_alternating_simulation(problem)
    heuristic = AS.build_heuristic(symmodel,[qT])
    AS.plot_heuristic!(heuristic)
    cells = Cell[]
    for (i,rec) in enumerate(rectangles)
        outneighbors = LightGraphs.outneighbors(symmodel.autom,i)
        n = length(outneighbors)
        coarse_abstraction =  nothing
        fine_abstraction = nothing
        heuristics = Dict{Int, Any}()
        controllers = Dict{NTuple{2,Int}, Any}()
        lower_bound = heuristic.dists[i]
        upper_bound = Inf
        local_target_set = Dict{Int, Any}()
        local_init_set = Dict{Int, Any}()
        cell = Cell(i,rec,outneighbors,coarse_abstraction,fine_abstraction,heuristics,controllers,
                    lower_bound,upper_bound,local_target_set,local_init_set)
        push!(cells,cell)
    end
    return cells
end

function compute_local_sets!(cells,contsys,Udom,q0,_I_,qT,_T_,compute_reachable_set)
    cells[q0].local_init_set[-1] = _I_
    cells[qT].local_target_set[-1] = _T_
    for cell in cells
        R = compute_reachable_set(cell.hyperrectangle,contsys,Udom)
        for i in cell.outneighbors
            neigh = cells[i]
            I1 = AB.intersect(R, neigh.hyperrectangle)
            neigh.local_init_set[cell.index] = I1
            cell.local_target_set[neigh.index] = I1
        end
    end
end

function plot_local_set!(cell)
    for neigh_i=-1:11
        if getkey(cell.local_target_set,neigh_i,false) != false
            println(neigh_i)
            h = cell.local_target_set[neigh_i]
            plot!(U.rectangle(h.lb,h.ub), opacity=.4,color=:blue)
        end
    end
end

end # end module
