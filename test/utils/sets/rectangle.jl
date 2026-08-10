module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Random
import LazySets
using StaticArrays: SVector

@testset "Hyperrectangle from bounds" begin
    R = LazySets.Hyperrectangle(; low = [1.0, 2.0], high = [3.0, 4.0])
    @test all(LazySets.low(R) .== [1.0, 2.0])
    @test all(LazySets.high(R) .== [3.0, 4.0])

    # Why the repo always writes the keyword form: positionally, the arguments are the *centre
    # and the radius*, so the same numbers build a different box and nothing complains.
    Rpos = LazySets.Hyperrectangle([1.0, 2.0], [3.0, 4.0])
    @test all(LazySets.low(Rpos) .== [-2.0, -2.0])
    @test all(LazySets.high(Rpos) .== [4.0, 6.0])

    # SVector bounds give SVector storage — the internal `_Box{N,T}` containers
    # (`set_in_period`, the certifier's tube) are typed on exactly this.
    Rs = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(1.0, 1.0))
    @test Rs isa UT._Box{2, Float64}

    # crossed bounds throw (an empty region is LazySets.EmptySet, never a box)
    @test_throws Exception LazySets.Hyperrectangle(; low = [1.0, 1.0], high = [-1.0, -1.0])
end

@testset "Ecosystem verbs (LazySets)" begin
    R1 = LazySets.Hyperrectangle(; low = [1.0, 1.0], high = [2.0, 2.0])
    R2 = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [3.0, 3.0])
    R5 = LazySets.Hyperrectangle(; low = [3.0, 3.0], high = [4.0, 4.0])

    @test R1 ⊆ R2
    @test !(R2 ⊆ R5)
    @test !isdisjoint(R1, R2)
    @test isdisjoint(R1, R5)

    @test LazySets.volume(R1) == 1
    @test LazySets.volume(R2) == 9

    S = LazySets.scale(0.5, R1)
    @test all(LazySets.low(S) .== [0.5, 0.5])
    @test all(LazySets.high(S) .== [1.0, 1.0])
end

@testset "Membership (points and boundaries)" begin
    R = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [1.0, 2.0])

    @test [0.0, 0.0] ∈ R
    @test [1.0, 2.0] ∈ R  # boundary
    @test [1.0, 2.0000001] ∉ R
    # LazySets membership is tolerance-based (rtol ≈ 1.5e-8), so use a point
    # clearly outside the set.
    @test [-1e-6, 1.0] ∉ R

    Rin = LazySets.Hyperrectangle(; low = [0.2, 0.3], high = [0.8, 1.9])
    @test Rin ⊆ R
    @test !(R ⊆ Rin)
end

@testset "Concrete intersection returns box or EmptySet" begin
    A = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [1.0, 1.0])

    # overlap
    B = LazySets.Hyperrectangle(; low = [0.5, 0.5], high = [2.0, 2.0])
    I = LazySets.intersection(A, B)
    @test all(LazySets.low(I) .== [0.5, 0.5])
    @test all(LazySets.high(I) .== [1.0, 1.0])
    @test !isempty(I)

    # disjoint -> EmptySet
    C = LazySets.Hyperrectangle(; low = [2.0, 2.0], high = [3.0, 3.0])
    I2 = LazySets.intersection(A, C)
    @test I2 isa LazySets.EmptySet
    @test isempty(I2)

    # touching at edge -> zero volume but not empty
    D = LazySets.Hyperrectangle(; low = [1.0, 0.0], high = [2.0, 1.0])
    I3 = LazySets.intersection(A, D)
    @test !isempty(I3)
    @test LazySets.volume(I3) == 0.0
end

@testset "Geometry helpers (center/h/r/dims/volume consistency)" begin
    R = LazySets.Hyperrectangle(; low = [1.0, 2.0], high = [3.0, 6.0])
    @test all(LazySets.center(R) .== [2.0, 4.0])
    @test all(UT.get_h(R) .== [2.0, 4.0])
    @test all(LazySets.radius_hyperrectangle(R) .== [1.0, 2.0])
    @test LazySets.dim(R) == 2

    @test LazySets.volume(R) == 8.0
end

@testset "Vertices" begin
    R = LazySets.Hyperrectangle(; low = [0.0, 10.0], high = [2.0, 12.0])
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
    R = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [1.0, 2.0])

    x = UT.sample(R)
    @test length(x) == 2
    @test x ∈ R

    xs = UT.samples(R, 200)
    @test length(xs) == 200
    @test all(z -> (z ∈ R), xs)
end

@testset "Periodic split (set_in_period)" begin
    rect_nowrap = LazySets.Hyperrectangle(; low = [0.2, 0.0], high = [0.8, 1.0])
    periodic_dims = SVector(1)     # dim 1 periodic
    periods = SVector(1.0)
    start = SVector(0.0)

    U1 = UT.set_in_period(rect_nowrap, periodic_dims, periods, start)
    @test length(U1.array) == 1
    @test U1.array[1] == rect_nowrap

    # wrapping case: [0.8, 1.2] wraps over period 1.0 -> two intervals
    rect_wrap = LazySets.Hyperrectangle(; low = [0.8, 0.0], high = [1.2, 1.0])
    U2 = UT.set_in_period(rect_wrap, periodic_dims, periods, start)
    @test length(U2.array) == 2

    # covers full period: width >= period -> single [start, start+period]
    rect_full = LazySets.Hyperrectangle(; low = [-0.5, 0.0], high = [1.5, 1.0]) # width 2.0 >= 1.0
    U3 = UT.set_in_period(rect_full, periodic_dims, periods, start)
    @test length(U3.array) == 1
    @test all(LazySets.low(U3.array[1]) .== [0.0, 0.0])
    @test all(LazySets.high(U3.array[1]) .== [1.0, 1.0])
end

end
