module TestMain

using Test
using Dionysos
using LinearAlgebra
using IntervalArithmetic
const DI = Dionysos
const UT = DI.Utils

function distance(E1, E2)
    return UT.pointCenterDistance(E1, E2.c)
end
function get_action(E1, E2)
    return (1.0, 1.0)
end

@testset "Tree" begin
    n_x = 2
    Ellispoids = [
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 8.0, [-10.0; -10.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 5.0, [0.0; -10.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 1.0, [-10.0; 0.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 3.0, [20.0; -10.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 3.0, [-1.0; 0.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 3.0, [1.0; -8.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 3.0, [-1.0; 5.0]),
        UT.Ellipsoid(Matrix{Float64}(I(n_x)) * 3.0, [3.0; 0.0]),
    ]
    tree = UT.Tree(Ellispoids[1])

    @test UT.get_nLeaves(tree) == 1
    @test UT.get_nNodes(tree) == 1
    @test UT.is_leave(tree.root) == true

    NclosestNodes, dists = UT.kNearestNeighbors(tree, Ellispoids[2], distance; k = 1)

    @test length(NclosestNodes) == 1
    @test dists[1] == 10.0

    action, cost = get_action(Ellispoids[2], NclosestNodes[1].state)
    node1 = UT.add_node!(tree, Ellispoids[2], NclosestNodes[1], action, cost)
    @test UT.is_leave(tree.root) == false
    @test UT.is_leave(node1) == true

    @test UT.get_nLeaves(tree) == 1
    @test UT.get_nNodes(tree) == 2

    NclosestNodes, dists = UT.kNearestNeighbors(tree, Ellispoids[2], distance; k = 2)
    @test length(NclosestNodes) == 2
    @test dists[1] == 0.0
    @test dists[2] == 10.0

    nNode3 = UT.add_closest_node!(tree, Ellispoids[3], distance, get_action)
    nNode4 = UT.add_closest_node!(tree, Ellispoids[4], distance, get_action)
    nNode5 = UT.add_closest_node!(tree, Ellispoids[5], distance, get_action)
    nNode6 = UT.add_closest_node!(tree, Ellispoids[6], distance, get_action)
    nNode7 = UT.add_closest_node!(tree, Ellispoids[7], distance, get_action)
    nNode8 = UT.add_closest_node!(tree, Ellispoids[8], distance, get_action)

    @test nNode3.path_cost == 1.0
    @test nNode5.path_cost == 2.0
    @test nNode7.path_cost == 3.0
    @test nNode8.path_cost == 3.0
    nNode3.path_cost = 3.0
    UT.propagate_cost_to_leaves(nNode3)
    @test nNode3.path_cost == 3.0
    @test nNode5.path_cost == 4.0
    @test nNode7.path_cost == 5.0
    @test nNode8.path_cost == 5.0
    UT.rewire(tree, nNode5, nNode6, 1.0, 1.0)
    @test nNode5.path_cost == 3.0
    @test nNode7.path_cost == 4.0
    @test nNode8.path_cost == 4.0
    @test UT.is_leave(nNode3) == true
    @test UT.is_leave(nNode6) == false
end

println("End test")
end
