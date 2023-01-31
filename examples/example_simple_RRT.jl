include("../src/Dionysos.jl")
using .Dionysos
UT = Dionysos.Utils
using Plots, Colors
using LinearAlgebra

function distance(E1,E2)
    return UT.pointCenterDistance(E1, E2.c)
end
function get_action(E1,E2)
    return 1.0
end

n_x = 2
Ellispoids = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*8.0, [-10.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*5.0, [0.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, [-10.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [20.0;-10.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [0.0;0.0]),
              UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [1.0;-8.0])
              ]


p = plot(aspect_ratio=:equal)

rrt = UT.RRT(Ellispoids[1])

NclosestNodes, dists = UT.findNClosestNode(rrt, Ellispoids[2], distance;N=1)
UT.add_node!(rrt,Ellispoids[2],NclosestNodes[1],get_action(Ellispoids[2],NclosestNodes[1].state),dists[1])

UT.add_closest_node!(rrt, Ellispoids[3], distance, get_action)
UT.add_closest_node!(rrt, Ellispoids[4], distance, get_action)
UT.add_closest_node!(rrt, Ellispoids[5], distance, get_action)
UT.add_closest_node!(rrt, Ellispoids[6], distance, get_action)

UT.plot_RRT!(rrt)
# UT.plot_path!(rrt, rrt.treeLeaves[1])


display(p)
