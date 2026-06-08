using DataStructures
using Dionysos

import Dionysos.System as ST

"""
    dijkstra_constrained_backward(autom::IndexedAutomatonList, initial_nodes, final_nodes, transition_cost, Delta)

Finds the optimal path by searching backwards from the final destinations to the origins.
Enforces that consecutive edge weights differ by <= Delta, and that the first/last edges 
themselves have a weight <= Delta.
"""
function dijkstra_constrained_backward(
    autom::ST.IndexedAutomatonList,
    initial_nodes,
    initial_u,
    final_nodes,
    final_u,
    transition_cost::Dict{Int, <:Real},
    Delta::Real
)

    initial_set = isa(initial_nodes, Set) ? initial_nodes : Set(initial_nodes)
    final_set = isa(final_nodes, Set) ? final_nodes : Set(final_nodes)

    # Trivial case: If an initial node is already a final node, the cost is 0 (0-edge path)
    common_nodes = intersect(initial_set, final_set)
    if !isempty(common_nodes)
        return [first(common_nodes)], Int[], 0.0
    end

    # Priority Queue stores: state (q, u) => accumulated_weight
    # (q, u) means: we are at node `q`, and taking forward symbol `u` moves us to the next node.
    pq = PriorityQueue{Tuple{Int, Int}, Float64}()
    dist = Dict{Tuple{Int, Int}, Float64}()
    
    # Successor map for forward path reconstruction: state (q, u) => next_state (q_next, u_next)
    successor = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()

    # 1. Initialization (Boundary condition for the LAST edge of the path)
    # The last edge 'u' entering a final node must satisfy: transition_cost[u] <= Delta
    for final_q in final_set
        for (q, u) in ST.pre(autom, final_q)
            edge_weight = transition_cost[u]
            if abs(edge_weight - final_u) <= Delta
                prev_state = (q, u)
                if !haskey(dist, prev_state) || edge_weight < dist[prev_state]
                    dist[prev_state] = edge_weight
                    pq[prev_state] = edge_weight
                    successor[prev_state] = (final_q, final_u)
                end
            end
        end
    end

    optimal_start_state = nothing
    min_total_weight = Inf

    # 2. Main Dijkstra Loop
    while !isempty(pq)
        curr_state, curr_dist = dequeue_pair!(pq)
        next_q, next_u = curr_state

        if curr_dist > dist[curr_state]
            continue
        end

        # 3. Termination Condition (Boundary condition for the FIRST edge of the path)
        # If q' is an initial node, and 'next_u' is close enough to 'initial_u', then it is the first edge of the forward path.
        if next_q in initial_set && abs(transition_cost[next_u] - transition_cost[initial_u]) <= Delta
            optimal_start_state = curr_state
            min_total_weight = curr_dist
            break
        end

        # 4. Backward Search
        for (q, u) in ST.pre(autom, next_q)
            edge_weight = transition_cost[u]
            if abs(edge_weight - transition_cost[next_u]) <= Delta
                next_dist = curr_dist + edge_weight
                prev_state = (q, u)
                if !haskey(dist, prev_state) || next_dist < dist[prev_state]
                    dist[prev_state] = next_dist
                    pq[prev_state] = next_dist
                    successor[prev_state] = curr_state
                end
            end
        end
    end

    if optimal_start_state === nothing
        return nothing, nothing, Inf
    end

    # 5. Path Reconstruction
    path_nodes = Int[]
    path_symbols = Int[]
    curr = optimal_start_state

    push!(path_nodes, curr[1])
    while curr[2] != 0
        push!(path_symbols, curr[2])
        curr = successor[curr]
        push!(path_nodes, curr[1])
    end

    return path_nodes, path_symbols, min_total_weight
end