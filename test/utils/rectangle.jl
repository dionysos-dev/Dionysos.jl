module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "RectangleOperator" begin
    R1 = UT.HyperRectangle([1,1],[2,2])
    R2 = UT.HyperRectangle([0,0],[3,3])
    X = R1 ∈ R2  
    @test X == true
    @test isempty(R1) == false
    R3 = UT.HyperRectangle([1,1],[-1,-1])
    @test isempty(R3) == true
    @test UT.volume(R1) == 1
    @test UT.volume(R2) == 9
    @test UT.volume(R3) == 0
    R4 = UT.HyperRectangle([1,1],[1,1])
    @test UT.volume(R4) == 0

    @test UT.is_intersection(R1,R2) == true
    R5 = UT.HyperRectangle([3,3],[4,4])
    @test UT.is_intersection(R1,R5) == false

    @test issubset(R1,R2) == true
    @test issubset(R2,R5) == false

    @test UT.scale(R1,0.5) == UT.HyperRectangle([0.5,0.5],[1.0,1])





end

sleep(0.1) # used for good printing
println("End test")
end