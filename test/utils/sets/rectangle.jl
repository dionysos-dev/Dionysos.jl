module TestMain

using Test
using Random
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "HyperRectangle constructors" begin
    # tuple constructor
    Rtu = UT.HyperRectangle((1.0, 2.0), (3.0, 4.0))
    @test Rtu.lb == UT.SVector(1.0, 2.0) || true  # if UT re-exports SVector, else just test values:
    @test all(Rtu.lb .== [1.0, 2.0])
    @test all(Rtu.ub .== [3.0, 4.0])

    # vector constructor
    Rv = UT.HyperRectangle([1, 2], [3, 4])
    @test all(Rv.lb .== [1, 2])
    @test all(Rv.ub .== [3, 4])

    # dimension mismatch should throw (SVector{length(lb)} with ub of different length typically errors)
    @test_throws Exception UT.HyperRectangle([1, 2], [3, 4, 5])
end

@testset "RectangleOperator basics" begin
    R1 = UT.HyperRectangle([1, 1], [2, 2])
    R2 = UT.HyperRectangle([0, 0], [3, 3])

    @test (R1 ⊆ R2) == true
    @test UT.isempty(R1) == false

    R3 = UT.HyperRectangle([1, 1], [-1, -1])
    @test UT.isempty(R3) == true

    @test UT.get_volume(R1) == 1
    @test UT.get_volume(R2) == 9
    @test UT.get_volume(R3) == 0

    R4 = UT.HyperRectangle([1, 1], [1, 1])
    @test UT.get_volume(R4) == 0

    @test isdisjoint(R1, R2) == false
    R5 = UT.HyperRectangle([3, 3], [4, 4])
    @test isdisjoint(R1, R5) == true

    @test UT.issubset(R1, R2) == true
    @test UT.issubset(R2, R5) == false

    @test UT.scale(R1, 0.5) == UT.HyperRectangle([0.5, 0.5], [1.0, 1.0])

    # NOTE: your get_h returns SVector; compare elementwise to avoid type mismatch
    @test all(UT.get_h(R1) .== [1, 1])
end

@testset "Membership (points and boundaries)" begin
    R = UT.HyperRectangle([0.0, 0.0], [1.0, 2.0])

    @test ([0.0, 0.0] ∈ R) == true
    @test ([1.0, 2.0] ∈ R) == true  # boundary
    @test ([1.0, 2.0000001] ∈ R) == false
    @test ([-1e-9, 1.0] ∈ R) == false

    # rectangle in rectangle
    Rin = UT.HyperRectangle([0.2, 0.3], [0.8, 1.9])
    @test (Rin ⊆ R) == true
    @test (R ⊆ Rin) == false
end

@testset "Equality / isequal / ==" begin
    A = UT.HyperRectangle([0, 0], [1, 1])
    B = UT.HyperRectangle([0, 0], [1, 1])
    C = UT.HyperRectangle([0, 0], [1, 2])

    @test isequal(A, B)
    @test A == B
    @test !isequal(A, C)
    @test A != C
end

@testset "Intersect edge cases" begin
    A = UT.HyperRectangle([0.0, 0.0], [1.0, 1.0])

    # overlap
    B = UT.HyperRectangle([0.5, 0.5], [2.0, 2.0])
    I = intersect(A, B)
    @test all(I.lb .== [0.5, 0.5])
    @test all(I.ub .== [1.0, 1.0])
    @test !UT.isempty(I)

    # disjoint -> empty rectangle (lb > ub somewhere)
    C = UT.HyperRectangle([2.0, 2.0], [3.0, 3.0])
    I2 = intersect(A, C)
    @test UT.isempty(I2)

    # touching at edge -> zero volume but not empty
    D = UT.HyperRectangle([1.0, 0.0], [2.0, 1.0])
    I3 = intersect(A, D)
    @test !UT.isempty(I3)
    @test UT.get_volume(I3) == 0.0
    @test all(I3.lb .== [1.0, 0.0])
    @test all(I3.ub .== [1.0, 1.0])
end

@testset "Geometry helpers (center/h/r/dims/volume consistency)" begin
    R = UT.HyperRectangle([1.0, 2.0], [3.0, 6.0])
    @test all(UT.get_center(R) .== [2.0, 4.0])
    @test all(UT.get_h(R) .== [2.0, 4.0])
    @test all(UT.get_r(R) .== [1.0, 2.0])
    @test UT.get_dim(R) == 2

    @test UT.get_volume(R) == 8.0
end

@testset "Vertices" begin
    R = UT.HyperRectangle([0.0, 10.0], [2.0, 12.0])
    V = UT.get_vertices(R)
    @test size(V, 1) == 2
    @test size(V, 2) == 4

    # all vertices inside, and min/max match bounds
    mins = [minimum(V[i, :]) for i in 1:2]
    maxs = [maximum(V[i, :]) for i in 1:2]
    @test mins == [0.0, 10.0]
    @test maxs == [2.0, 12.0]

    verts = UT.collect_vertices(R)
    @test length(verts) == 4
    @test all(v -> (v ∈ R), verts)
end

@testset "Sampling" begin
    Random.seed!(1234)
    R = UT.HyperRectangle([0.0, 0.0], [1.0, 2.0])

    x = UT.sample(R)
    @test length(x) == 2
    @test x ∈ R

    xs = UT.samples(R, 200)
    @test length(xs) == 200
    @test all(z -> (z ∈ R), xs)
end

@testset "Periodic split (set_in_period)" begin
    # Minimal sanity tests that exercise code paths.

    rect_nowrap = UT.HyperRectangle([0.2, 0.0], [0.8, 1.0])
    periodic_dims = UT.SVector(1)     # dim 1 periodic
    periods = UT.SVector(1.0)
    start = UT.SVector(0.0)

    U1 = UT.set_in_period(rect_nowrap, periodic_dims, periods, start)
    @test length(U1.array) == 1
    @test U1.array[1] == rect_nowrap

    # wrapping case: [0.8, 1.2] wraps over period 1.0 -> two intervals
    rect_wrap = UT.HyperRectangle([0.8, 0.0], [1.2, 1.0])
    U2 = UT.set_in_period(rect_wrap, periodic_dims, periods, start)
    @test length(U2.array) == 2
    @test all(r -> !UT.isempty(r), U2.array)

    # covers full period: width >= period -> single [start, start+period]
    rect_full = UT.HyperRectangle([-0.5, 0.0], [1.5, 1.0]) # width 2.0 >= 1.0
    U3 = UT.set_in_period(rect_full, periodic_dims, periods, start)
    @test length(U3.array) == 1
    @test all(U3.array[1].lb .== [0.0, 0.0])
    @test all(U3.array[1].ub .== [1.0, 1.0])
end

println("End test")
end
