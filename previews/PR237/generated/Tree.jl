using LinearAlgebra, Plots, Colors

using Dionysos
UT = Dionysos.Utils

function distance(E1, E2)
    return UT.pointCenterDistance(E1, E2.c)
end

function get_action(E1, E2)
    return (1.0, 1.0)
end

Ellipsoids = [UT.Ellipsoid(Matrix{Float64}(I(2))*8.0, [-10.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*5.0, [0.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*1.0, [-10.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*3.0, [20.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*3.0, [-1.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*3.0, [1.0;-8.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*3.0, [-1.0;5.0]),
              UT.Ellipsoid(Matrix{Float64}(I(2))*3.0, [3.0;0.0])]

tree = UT.Tree(Ellipsoids[1])

action, cost = get_action(Ellipsoids[2], tree.root.state)
nNode2 = UT.add_node!(tree, Ellipsoids[2], tree.root, action, cost)

nNode3 = UT.add_closest_node!(tree, Ellipsoids[3], distance, get_action)

nNode4 = UT.add_closest_node!(tree, Ellipsoids[4], distance, get_action)
nNode5 = UT.add_closest_node!(tree, Ellipsoids[5], distance, get_action)
nNode6 = UT.add_closest_node!(tree, Ellipsoids[6], distance, get_action)
nNode7 = UT.add_closest_node!(tree, Ellipsoids[7], distance, get_action)
nNode8 = UT.add_closest_node!(tree, Ellipsoids[8], distance, get_action)

UT.print_data(tree)
fig = plot(aspect_ratio=:equal)
plot!(tree; arrowsB=true, cost=true)

nNode3.path_cost = 5.0
UT.propagate_cost_to_leaves(nNode3)
fig = plot(aspect_ratio=:equal)
plot!(tree)

UT.rewire(tree, nNode5, nNode6, 1.0, 1.0)
fig = plot(aspect_ratio=:equal)
plot!(tree)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

