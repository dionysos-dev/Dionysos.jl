using Dionysos

import Dionysos.System as ST

include("constrained_Dijkstra.jl")

graph = ST.NewIndexedAutomatonList(9, 4)

# Add transitions:
#  - 1 = UP
#  - 2 = DOWN
#  - 3 = LEFT
#  - 4 = RIGHT

for u in 1:4
    for q in 4:9
        ST.add_transition!(graph, q, q-3, 1)
    end
    for q in 1:6
        ST.add_transition!(graph, q, q+3, 2)
    end
    for q in [2, 3, 5, 6, 8, 9]
        ST.add_transition!(graph, q, q-1, 3)
    end
    for q in [1, 2, 4, 5, 7, 8]
        ST.add_transition!(graph, q, q+1, 4) end
end

# Define input vector associated with transition:
#  - UP = [0, 1]
#  - DOWN = [0, -2] (suppose it requires more force to go down)
#  - LEFT = [-1 0]
#  - RIGHT = [1 0]
transition_weight = Dict([
    (0, 0.0),
    (1, 1.0),
    (2, 2.0),
    (3, 1.0),
    (4, 1.0)
])

dijkstra_constrained_backward(graph, 1, 0, 9, 0, transition_weight, 1)