include("../src/Dionysos.jl")
using .Dionysos
UT = Dionysos.Utils
using Plots, Colors
using LinearAlgebra

function distance(E1,E2)
    return UT.pointCenterDistance(E1, E2.c)
end
function get_action(E1,E2)
    return (1.0, 1.0)
end

n_x = 2
Ellispoids = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*8.0, [-10.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*5.0, [0.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, [-10.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [20.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [-1.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [1.0;-8.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [-1.0;5.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [3.0;0.0])
              ]


p = plot(aspect_ratio=:equal)

tree = UT.Tree(Ellispoids[1])
println(UT.is_leave(tree.root))

NclosestNodes, dists = UT.kNearestNeighbors(tree, Ellispoids[2], distance; k=1)

action, cost = get_action(Ellispoids[2],NclosestNodes[1].state)
nNode1 = UT.add_node!(tree,Ellispoids[2],NclosestNodes[1], action, cost)
println(UT.is_leave(tree.root))

nNode2 = UT.add_closest_node!(tree, Ellispoids[3], distance, get_action)
nNode3 = UT.add_closest_node!(tree, Ellispoids[4], distance, get_action)
nNode4 = UT.add_closest_node!(tree, Ellispoids[5], distance, get_action)
nNode5 = UT.add_closest_node!(tree, Ellispoids[6], distance, get_action)
nNode6 = UT.add_closest_node!(tree, Ellispoids[7], distance, get_action)
nNode7 = UT.add_closest_node!(tree, Ellispoids[8], distance, get_action)

println(nNode2.path_cost)
nNode2.path_cost = 3.0
UT.propagate_cost_to_leaves(nNode2)

UT.rewire(tree, nNode4, nNode5, 1.0, 1.0)
UT.plot_Tree!(tree)

# UT.plot_path!(tree.leaves[1])
UT.print_data(tree)


display(p)
