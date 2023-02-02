module TestMain

using Test
using Dionysos
using LinearAlgebra
using IntervalArithmetic
const DI = Dionysos
const UT = DI.Utils

function distance(E1,E2)
    return UT.pointCenterDistance(E1, E2.c)
end
function get_action(E1,E2)
    return 1.0
end

@testset "Tree" begin
    n_x = 2
    Ellispoids = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*8.0, [-10.0;-10.0]),
                UT.Ellipsoid(Matrix{Float64}(I(n_x))*5.0, [0.0;-10.0]),
                UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, [-10.0;0.0]),
                UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [20.0;-10.0]),
                UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [0.0;0.0]),
                UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [1.0;-8.0])
                ]
    rrt = UT.RRT(Ellispoids[1])

    @test UT.get_nLeaves(rrt) == 1
    @test UT.get_nNodes(rrt) == 1

    NclosestNodes, dists = UT.findNClosestNode(rrt, Ellispoids[2], distance;N=1)

    @test length(NclosestNodes) == 1
    @test dists[1] == 10.0

    UT.add_node!(rrt,Ellispoids[2],NclosestNodes[1],get_action(Ellispoids[2],NclosestNodes[1].state),dists[1])

    @test UT.get_nLeaves(rrt) == 1
    @test UT.get_nNodes(rrt) == 2

    NclosestNodes, dists = UT.findNClosestNode(rrt, Ellispoids[2], distance;N=2)
    @test length(NclosestNodes) == 2
    @test dists[1] == 0.0
    @test dists[2] == 10.0
end


sleep(0.1) # used for good printing
println("End test")
end