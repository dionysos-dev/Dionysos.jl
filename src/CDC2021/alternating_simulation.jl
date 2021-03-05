module AlternatingSimulation
"""
    Build alternating simulation.
"""

using ..Abstraction
const AB = Abstraction

using ..Utils
const U = Utils

using ..DomainList
const D = DomainList

using StaticArrays, LightGraphs, SimpleWeightedGraphs, Plots

abstract type AbstractionProblem end

struct symmodelProblem{M,S,R,T,G,E} <: AbstractionProblem
    symmodel::M
    contsys::S
    compute_reachable_set::R     # function to compute an hyperrectangle overpproximation of the reachable set in 1 time step
    minimum_transition_cost::T   # function to compute the minimum transition cost between two cells in 1 time step
    get_possible_transitions::G  # function to compute a list of potential neighbour cells.
    ext::E
end

function symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,get_possible_transitions;ext=nothing)
    return symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,get_possible_transitions,ext)
end

function get_ncells(problem::symmodelProblem)
    symmodel = problem.symmodel
    return AB.get_ncells(symmodel.Xdom)
end
function enum_cells(problem::symmodelProblem)
    return [i for i=1:AB.get_ncells(problem.symmodel.Xdom)]
end

function build_alternating_simulation(problem::symmodelProblem)
    digraph = SimpleWeightedDiGraph(get_ncells(problem))
    for source in enum_cells(problem)
        for neighbor in problem.get_possible_transitions(problem,source)
            cost = problem.minimum_transition_cost(problem.symmodel,problem.contsys,source,neighbor)
            if cost < Inf
                add_edge!(digraph, source, neighbor, cost)
            end
        end
    end
    return digraph
end

# return the list of possible neighboring cells (in one time step) for the big cell abstraction
function get_possible_transitions_1(problem::symmodelProblem,source::Int)
    hyperrectangle = problem.symmodel.Xdom.elems[source]
    reachable_sets = problem.compute_reachable_set(hyperrectangle,problem.contsys,problem.symmodel.Udom)
    neighbors = Int[]
    for (index,neigh) in enumerate(problem.symmodel.Xdom.elems)
        if source!=index && U.is_intersection(reachable_sets,neigh)
            push!(neighbors,index)
        end
    end
    return neighbors
end
# return the list of possible neighboring cells (in one time step) for the coarse abstraction
function get_possible_transitions_2(problem::symmodelProblem,source::Int)
    symmodel = problem.symmodel
    pos = AB.get_xpos_by_state(symmodel,source)
    grid = symmodel.Xdom.grid
    c = AB.get_coord_by_pos(grid,pos)
    h = grid.h
    hyperrectangle = AB.HyperRectangle(c-h/2,c+h/2)
    reachable_set = problem.compute_reachable_set(hyperrectangle,problem.contsys,problem.symmodel.Udom)
    println("reachable set computed done")

    symbols = U.get_symbol(symmodel,reachable_set,AB.OUTER)
    return symbols
end=#

# return the list of possible neighboring cells (in one time step) for the coarse abstraction
function get_possible_transitions_2(problem::symmodelProblem,source::Int)
    symmodel = problem.symmodel
    Xdom = symmodel.Xdom
    grid = Xdom.grid
    pos = AB.get_xpos_by_state(symmodel,source)
    c = AB.get_coord_by_pos(grid,pos)
    h = grid.h
    hyperrectangle = AB.HyperRectangle(c-h/2,c+h/2)
    reachable_set = problem.compute_reachable_set(hyperrectangle,problem.contsys,problem.symmodel.Udom)
    # il se peut que pour les dimensions non periodiques, que les rectangles soient en dehors du domaine si
    # l'overapprox du reachable set est tres grande. (si c'est le cas, le nombre de cell à enumemer peut etre immense,
    # alors que la plupard sont hors du domaine) d'où lims dans general_domain
    reachable_sets = D.set_rec_in_period(Xdom.periodic,Xdom.periods,Xdom.T0,reachable_set)
    symbols = U.get_symbols(symmodel,reachable_sets,AB.OUTER)
    return symbols
end

function get_possible_transitions_3(problem::symmodelProblem,source::Int)
    symmodel = problem.symmodel
    Xdom = symmodel.Xdom
    reachable_set = problem.ext[source]
    # il se peut que pour les dimensions non periodiques, que les rectangles soient en dehors du domaine si
    # l'overapprox du reachable set est tres grande. (si c'est le cas, le nombre de cell à enumemer peut etre immense,
    # alors que la plupard sont hors du domaine) d'où lims dans general_domain
    reachable_sets = D.set_rec_in_period(Xdom.periodic,Xdom.periods,Xdom.T0,reachable_set)
    symbols = U.get_symbols(symmodel,reachable_sets,AB.OUTER)
    return symbols
end

## Heuristic
abstract type Heuristic end

struct symmodelHeuristic <: Heuristic
    symmodel
    dists
end

function build_heuristic(symmodel,initlist)
    result = LightGraphs.dijkstra_shortest_paths(symmodel.autom,initlist)
    heuristic = symmodelHeuristic(symmodel,result.dists)
    return heuristic
end

function get_min_value_heurisitic(heuristic,subset)
    symmodel = heuristic.symmodel
    Xdom = symmodel.Xdom
    grid = Xdom.grid
    subdomain = AB.DomainList(grid)
    AB.add_subset!(subdomain, Xdom, subset, AB.OUTER)
    val = Inf
    for pos in AB.enum_pos(subdomain)
        val = min(val,heuristic.dists[AB.get_state_by_xpos(symmodel, pos)])
    end
    return val
end

function get_coord(Xdom::AB.DomainList, pos)
    return AB.get_coord_by_pos(Xdom.grid,pos)
end
function get_coord(Xdom::U.CustomList, rec)
    return U.center(rec)
end

function plot_heuristic!(heuristic::symmodelHeuristic;dims=[1,2],opacity=0.2,color=:red)
    symmodel = heuristic.symmodel
    dists = heuristic.dists
    U.plot_domain!(symmodel.Xdom,dims=dims,opacity=opacity,color=color)
    i = 1
    for elem in AB.enum_pos(symmodel.Xdom)
        x = get_coord(symmodel.Xdom,elem)[dims]
        if dists[i] != Inf
            annotate!(x[1], x[2], text(Int(dists[i]),4), :color)
        else
            annotate!(x[1], x[2], text(Inf,4), :color)
        end
        i = i+1
    end
end

end # module
