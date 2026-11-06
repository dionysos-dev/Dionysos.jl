# # Tree
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Tree.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Tree.ipynb)
#
# In this file, we present the Tree data structure, which is a tree composed of costs on edges and an underlying metric between node states.
# In this simple example, the states are two-dimensional ellipsoids.
# 
# First, let us import a few packages that are necessary to run this example.
using LinearAlgebra, Plots, Colors

# The main package [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) provides most important data structures that we will need.
using Dionysos
const UT = Dionysos.Utils

# We define the underlying metric between node states
distance(E1::UT.Ellipsoid, E2::UT.Ellipsoid) = UT.pointCenterDistance(E1, E2.c)

# We define the action function to compute a transition between two states
get_action(E1::UT.Ellipsoid, E2::UT.Ellipsoid) = (1.0, 1.0)

# We define the ellipsoids that will make up our tree states
Ellipsoids = [
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 8.0, [-10.0; -10.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 5.0, [0.0; -10.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 1.0, [-10.0; 0.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 3.0, [20.0; -10.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 3.0, [-1.0; 0.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 3.0, [1.0; -8.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 3.0, [-1.0; 5.0]),
    UT.Ellipsoid(Matrix{Float64}(I(2)) * 3.0, [3.0; 0.0]),
]

# We define the root of the tree
tree = UT.Tree(Ellipsoids[1])

# Compute the transition between Ellipsoids[2] and the root of the tree
action, cost = get_action(Ellipsoids[2], tree.root.state)
nNode2 = UT.add_node!(tree, Ellipsoids[2], tree.root, action, cost)

# Connect Ellipsoids[3] to its closest node according to the underlying metric
nNode3 = UT.add_closest_node!(tree, Ellipsoids[3], distance, get_action)

# Connect the other ellipsoids
nNode4 = UT.add_closest_node!(tree, Ellipsoids[4], distance, get_action)
nNode5 = UT.add_closest_node!(tree, Ellipsoids[5], distance, get_action)
nNode6 = UT.add_closest_node!(tree, Ellipsoids[6], distance, get_action)
nNode7 = UT.add_closest_node!(tree, Ellipsoids[7], distance, get_action)
nNode8 = UT.add_closest_node!(tree, Ellipsoids[8], distance, get_action)

# Plot the tree
println(tree)
fig = plot(; aspect_ratio = :equal)
plot!(tree; arrowsB = true, cost = true)

# We change the node's cost and update the tree accordingly
nNode3.path_cost = 5.0
UT.propagate_cost_to_leaves(nNode3)
fig = plot(; aspect_ratio = :equal)
plot!(tree)

# We change the node's parent
UT.rewire(tree, nNode5, nNode6, 1.0, 1.0)
fig = plot(; aspect_ratio = :equal)
plot!(tree)
