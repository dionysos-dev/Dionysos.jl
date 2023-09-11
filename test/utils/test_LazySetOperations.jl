module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "LazySetOperations" begin
    A1 = UT.HyperRectangle([-1.0 -1.0], [1.0 1.0])
    c1 = [0.0; 0.0]
    P = [
        1.0 0.0
        0.0 1.0
    ]
    B1 = UT.Ellipsoid(P, c1)
    AminB1 = UT.LazySetMinus(A1, B1)
    B = UT.get_B(AminB1)
    @test UT.get_center(B) == c1
    @test UT.get_shape(B) == P

    D = UT.LazyUnionSetArray([])
    @test UT.isempty(D)
    @test UT.isempty(UT.get_sets(D))

    c2 = [1.0; 1.0]
    B2 = UT.Ellipsoid(P, c2)
    C = UT.LazyUnionSetArray([B1, B2])
    Csets = UT.get_sets(C)
    @test Csets[1] == B1
    @test Csets[2] == B2

    A2 = UT.HyperRectangle([1.0 1.0], [2.0 2.0])
    E = UT.LazyUnionSetArray([A1, A2])
    A3 = UT.HyperRectangle([0.5 0.5], [1.5 1.5])
    I = UT.get_sets(intersect(E, A3))
    @test I[1] == UT.HyperRectangle([0.5 0.5], [1.0 1.0])
    @test I[2] == UT.HyperRectangle([1.0 1.0], [1.5 1.5])
end

println("End test")
end
