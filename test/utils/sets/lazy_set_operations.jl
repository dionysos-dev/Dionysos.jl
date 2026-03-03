module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "LazySetOperations" begin
    # ----------------------------
    # Basic Minus
    # ----------------------------
    A1 = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])

    c1 = [0.0, 0.0]
    P = [1.0 0.0; 0.0 1.0]
    B1 = UT.Ellipsoid(P, c1)

    AminB1 = UT.LazySetMinus(A1, B1)

    @test AminB1.B == B1
    @test UT.get_center(AminB1.B) == c1
    @test UT.get_shape(AminB1.B) == P

    # ----------------------------
    # Empty union
    # ----------------------------
    U0 = UT.LazySetUnion{2, Float64}()
    @test isempty(U0)
    @test isempty(UT.get_sets(U0))

    # ----------------------------
    # Union stores nodes
    # ----------------------------
    c2 = [1.0, 1.0]
    B2 = UT.Ellipsoid(P, c2)

    Ue = UT.LazySetUnion([B1, B2])
    Ue_sets = UT.get_sets(Ue)

    @test Ue_sets[1] == B1
    @test Ue_sets[2] == B2

    # ----------------------------
    # Intersection distribution over union
    # ----------------------------
    A2 = UT.HyperRectangle([1.0, 1.0], [2.0, 2.0])
    Urect = UT.LazySetUnion([A1, A2])

    A3 = UT.HyperRectangle([0.5, 0.5], [1.5, 1.5])

    I = intersect(Urect, A3)
    Isets = UT.get_sets(I)

    @test Isets[1] == UT.HyperRectangle([0.5, 0.5], [1.0, 1.0])
    @test Isets[2] == UT.HyperRectangle([1.0, 1.0], [1.5, 1.5])
end

println("End test")

end # module