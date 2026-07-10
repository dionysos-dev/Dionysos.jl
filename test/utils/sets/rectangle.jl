module TestMain

using Test
using Random
import LazySets
using StaticArrays: SVector
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "box constructors" begin
    # tuple constructor
    Rtu = UT.box((1.0, 2.0), (3.0, 4.0))
    @test all(LazySets.low(Rtu) .== [1.0, 2.0])
    @test all(LazySets.high(Rtu) .== [3.0, 4.0])

    # vector constructor promotes Int
    Rv = UT.box([1, 2], [3, 4])
    @test all(LazySets.low(Rv) .== [1.0, 2.0])
    @test all(LazySets.high(Rv) .== [3.0, 4.0])

    # SVector constructor keeps static storage
    Rs = UT.box(SVector(0.0, 0.0), SVector(1.0, 1.0))
    @test Rs isa UT.Box{2, Float64}

    # dimension mismatch throws
    @test_throws DimensionMismatch UT.box([1, 2], [3, 4, 5])

    # crossed bounds throw (an empty region is LazySets.EmptySet, never a box)
    @test_throws Exception UT.box([1.0, 1.0], [-1.0, -1.0])
end

@testset "Ecosystem verbs (LazySets)" begin
    R1 = UT.box([1, 1], [2, 2])
    R2 = UT.box([0, 0], [3, 3])
    R5 = UT.box([3, 3], [4, 4])

    @test R1 ⊆ R2
    @test !(R2 ⊆ R5)
    @test !isdisjoint(R1, R2)
    @test isdisjoint(R1, R5)

    @test UT.get_volume(R1) == 1
    @test UT.get_volume(R2) == 9

    S = LazySets.scale(0.5, R1)
    @test all(LazySets.low(S) .== [0.5, 0.5])
    @test all(LazySets.high(S) .== [1.0, 1.0])
end

@testset "Membership (points and boundaries)" begin
    R = UT.box([0.0, 0.0], [1.0, 2.0])

    @test [0.0, 0.0] ∈ R
    @test [1.0, 2.0] ∈ R  # boundary
    @test [1.0, 2.0000001] ∉ R
    # LazySets membership is tolerance-based (rtol ≈ 1.5e-8), so use a point
    # clearly outside the set.
    @test [-1e-6, 1.0] ∉ R

    Rin = UT.box([0.2, 0.3], [0.8, 1.9])
    @test Rin ⊆ R
    @test !(R ⊆ Rin)
end

@testset "Concrete intersection returns box or EmptySet" begin
    A = UT.box([0.0, 0.0], [1.0, 1.0])

    # overlap
    B = UT.box([0.5, 0.5], [2.0, 2.0])
    I = LazySets.intersection(A, B)
    @test all(LazySets.low(I) .== [0.5, 0.5])
    @test all(LazySets.high(I) .== [1.0, 1.0])
    @test !isempty(I)

    # disjoint -> EmptySet
    C = UT.box([2.0, 2.0], [3.0, 3.0])
    I2 = LazySets.intersection(A, C)
    @test I2 isa LazySets.EmptySet
    @test isempty(I2)

    # touching at edge -> zero volume but not empty
    D = UT.box([1.0, 0.0], [2.0, 1.0])
    I3 = LazySets.intersection(A, D)
    @test !isempty(I3)
    @test UT.get_volume(I3) == 0.0
end

@testset "Geometry helpers (center/h/r/dims/volume consistency)" begin
    R = UT.box([1.0, 2.0], [3.0, 6.0])
    @test all(UT.get_center(R) .== [2.0, 4.0])
    @test all(UT.get_h(R) .== [2.0, 4.0])
    @test all(UT.get_r(R) .== [1.0, 2.0])
    @test UT.get_dim(R) == 2

    @test UT.get_volume(R) == 8.0
end

@testset "Vertices" begin
    R = UT.box([0.0, 10.0], [2.0, 12.0])
    verts = LazySets.vertices_list(R)
    @test length(verts) == 4
    @test all(v -> (v ∈ R), verts)

    mins = [minimum(v[i] for v in verts) for i in 1:2]
    maxs = [maximum(v[i] for v in verts) for i in 1:2]
    @test mins == [0.0, 10.0]
    @test maxs == [2.0, 12.0]
end

@testset "Sampling" begin
    Random.seed!(1234)
    R = UT.box([0.0, 0.0], [1.0, 2.0])

    x = UT.sample(R)
    @test length(x) == 2
    @test x ∈ R

    xs = UT.samples(R, 200)
    @test length(xs) == 200
    @test all(z -> (z ∈ R), xs)
end

@testset "Periodic split (set_in_period)" begin
    rect_nowrap = UT.box([0.2, 0.0], [0.8, 1.0])
    periodic_dims = SVector(1)     # dim 1 periodic
    periods = SVector(1.0)
    start = SVector(0.0)

    U1 = UT.set_in_period(rect_nowrap, periodic_dims, periods, start)
    @test length(U1.array) == 1
    @test U1.array[1] == rect_nowrap

    # wrapping case: [0.8, 1.2] wraps over period 1.0 -> two intervals
    rect_wrap = UT.box([0.8, 0.0], [1.2, 1.0])
    U2 = UT.set_in_period(rect_wrap, periodic_dims, periods, start)
    @test length(U2.array) == 2

    # covers full period: width >= period -> single [start, start+period]
    rect_full = UT.box([-0.5, 0.0], [1.5, 1.0]) # width 2.0 >= 1.0
    U3 = UT.set_in_period(rect_full, periodic_dims, periods, start)
    @test length(U3.array) == 1
    @test all(LazySets.low(U3.array[1]) .== [0.0, 0.0])
    @test all(LazySets.high(U3.array[1]) .== [1.0, 1.0])
end

println("End test")
end
