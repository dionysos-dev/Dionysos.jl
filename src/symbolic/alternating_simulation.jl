using Graphs, SimpleWeightedGraphs

function symmodelAS(Xdom, Udom, sys, minimum_transition_cost, get_possible_transitions)
    symmodel = NewSymbolicModelListList(Xdom, Udom)
    ncells = DO.get_ncells(Xdom)
    digraph = SimpleWeightedDiGraph(ncells)
    for source in enum_cells(symmodel)
        for neighbor in get_possible_transitions(symmodel, sys, source)
            cost = minimum_transition_cost(symmodel, sys, source, neighbor)
            if cost < Inf
                SimpleWeightedGraphs.add_edge!(digraph, source, neighbor, cost)
            end
        end
    end
    symmodel = with_automaton(symmodel, digraph)
    return symmodel
end

function get_transitions_1(symmodel, sys, source::Int, compute_reachable_set)
    Xdom = symmodel.Xdom
    grid = Xdom.grid
    pos = get_xpos_by_state(symmodel, source)
    c = DO.get_coord_by_pos(grid, pos)
    h = grid.h
    hyperrectangle = UT.HyperRectangle(c - h / 2, c + h / 2)
    reachable_set = compute_reachable_set(hyperrectangle, sys, symmodel.Udom)
    # il se peut que pour les dimensions non periodiques, que les rectangles soient en dehors du domaine si
    # l'overapprox du reachable set est tres grande. (si c'est le cas, le nombre de cell à enumemer peut etre immense,
    # alors que la plupard sont hors du domaine) d'où lims dans general_domain
    reachable_sets =
        DO.set_rec_in_period(Xdom.periodic, Xdom.periods, Xdom.T0, reachable_set)
    symbols = get_symbols(symmodel, reachable_sets, DO.OUTER)
    return symbols
end

## Heuristic
abstract type Heuristic end

struct symmodelHeuristic <: Heuristic
    symmodel::Any
    dists::Any
end

function build_heuristic(symmodel, initlist)
    result = Graphs.dijkstra_shortest_paths(symmodel.autom, initlist)
    heuristic = symmodelHeuristic(symmodel, result.dists)
    return heuristic
end
