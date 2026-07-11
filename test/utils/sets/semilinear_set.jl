
module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using LazySets
using Polyhedra
using CDDLib
import LinearAlgebra as LA

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

const Poly = LazySets.HPolytope
const HS = LazySets.HalfSpace

backend = CDDLib.Library()

"""
Axis-aligned box in H-representation:
    lb[i] <= x[i] <= ub[i]
"""
function hbox(lb::AbstractVector, ub::AbstractVector)
    n = length(lb)
    cons = HS[]
    for i in 1:n
        ei = zeros(eltype(float(lb[i])), n)
        ei[i] = 1.0
        push!(cons, HS(copy(ei), ub[i]))     #  x_i <= ub_i
        push!(cons, HS(-copy(ei), -lb[i]))   # -x_i <= -lb_i  <=> x_i >= lb_i
    end
    return Poly(cons)
end

@testset "Poly helpers" begin
    P = hbox([0.0, 0.0], [1.0, 2.0])

    Q = UT._as_hpolytope(P)
    @test Q isa Poly

    B = LazySets.Hyperrectangle([0.5, 1.0], [0.5, 1.0])
    @test UT._as_hpolytope(B) isa Poly

    C = UT.clean_poly(P)
    @test C isa Poly

    @test LazySets.dim(P) == 2
    @test !isempty(P)

    E = LazySets.EmptySet{Float64}(2)
    @test isempty(E)
end

@testset "Poly intersection (LazySets-native)" begin
    P = hbox([0.0, 0.0], [1.0, 1.0])
    Q = hbox([0.5, 0.5], [2.0, 2.0])

    I = LazySets.intersection(P, Q)
    @test !isempty(I)

    @test [0.75, 0.75] ∈ I
    @test [0.25, 0.25] ∉ I
    @test [1.5, 1.5] ∉ I

    # disjoint polytopes: one LP, no intersection construction needed
    R = hbox([5.0, 5.0], [6.0, 6.0])
    @test UT.is_disjoint(P, R)
    @test !UT.is_disjoint(P, Q)
end

@testset "Poly volume (LazySets-native)" begin
    P = hbox([0.0, 0.0], [1.0, 2.0])
    v = LazySets.volume(P; backend = backend)
    @test isapprox(v, 2.0; atol = 1e-8)

    # in 2-D LazySets computes the area without a backend
    @test isapprox(LazySets.volume(P), 2.0; atol = 1e-8)
end

@testset "SemiLinearSet constructors and basics" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])

    S0 = UT.semilinear_set()
    @test S0 isa UT.SemiLinearSet
    @test isempty(S0.array)
    @test isempty(S0)

    S1 = UT.semilinear_set(P1)
    @test S1 isa UT.SemiLinearSet
    @test length(S1.array) == 1
    @test !isempty(S1)
    @test LazySets.dim(S1) == 2

    S2 = UT.semilinear_set([P1, P2])
    @test S2 isa UT.SemiLinearSet
    @test length(S2.array) == 2
    @test !isempty(S2)

    @test [0.5, 0.5] ∈ S2
    @test [2.5, 2.5] ∈ S2
    @test [1.5, 1.5] ∉ S2

    # non-polytopic parts are converted; a union of boxes is not a SemiLinearSet
    B = LazySets.Hyperrectangle([0.5, 0.5], [0.5, 0.5])
    @test UT.semilinear_set([B]) isa UT.SemiLinearSet
    @test !(LazySets.UnionSetArray([B]) isa UT.SemiLinearSet)
end

@testset "normalize_semilinear" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])

    # Empty polytope in H-representation: x <= 0 and x >= 1 simultaneously
    E = Poly([HS([1.0, 0.0], 0.0), HS([-1.0, 0.0], -1.0)])

    S = UT.semilinear_set([P1, E])
    Sn = UT.normalize_semilinear(S)

    @test length(Sn.array) == 1
    @test [0.5, 0.5] ∈ Sn
    @test [2.0, 2.0] ∉ Sn
end

@testset "SemiLinearSet nonemptiness" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])

    E = Poly([HS([1.0, 0.0], 0.0), HS([-1.0, 0.0], -1.0)])

    S1 = UT.semilinear_set([E])
    S2 = UT.semilinear_set([E, P1])

    @test isempty(S1)
    @test !isempty(S2)
end

@testset "SemiLinearSet intersection (LazySets-native)" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])
    Q = hbox([0.5, 0.5], [2.5, 2.5])

    S = UT.semilinear_set([P1, P2])

    I1 = LazySets.intersection(S, UT.semilinear_set(Q))
    @test I1 isa LazySets.UnionSetArray
    @test length(I1.array) == 2
    @test [0.75, 0.75] ∈ I1
    @test [2.25, 2.25] ∈ I1
    @test [1.5, 1.5] ∉ I1

    I2 = LazySets.intersection(S, Q)
    @test [0.75, 0.75] ∈ I2
    @test [2.25, 2.25] ∈ I2
    @test [1.5, 1.5] ∉ I2

    # empty parts are filtered: a set touching only one part gives that part
    Qlow = hbox([0.25, 0.25], [0.75, 0.75])
    I3 = LazySets.intersection(S, Qlow)
    @test !(I3 isa LazySets.UnionSetArray)  # single surviving part
    @test [0.5, 0.5] ∈ I3

    # fully disjoint -> EmptySet
    Qout = hbox([5.0, 5.0], [6.0, 6.0])
    @test LazySets.intersection(S, Qout) isa LazySets.EmptySet
end

@testset "poly_intersection" begin
    P = hbox([0.0, 0.0], [2.0, 2.0])
    Q = hbox([1.0, 1.0], [3.0, 3.0])
    I = UT.poly_intersection(P, Q)
    @test I isa Poly
    @test [1.5, 1.5] ∈ I

    R = hbox([5.0, 5.0], [6.0, 6.0])
    @test UT.poly_intersection(P, R) isa LazySets.EmptySet

    # touching at a face: degenerate intersection, must not throw
    T = hbox([2.0, 0.0], [3.0, 2.0])
    @test UT.poly_intersection(P, T) isa Union{Poly, LazySets.EmptySet}

    S = UT.semilinear_set([P, R])
    parts = UT.poly_intersection_parts(S, Q)
    @test length(parts) == 1
    @test [1.5, 1.5] ∈ parts[1]
end

@testset "Poly difference decomposition" begin
    P = hbox([0.0, 0.0], [2.0, 2.0])
    Q = hbox([0.5, 0.5], [1.5, 1.5])

    pieces = UT.set_difference_decompose(P, Q)
    @test pieces isa Vector{<:Poly}
    @test length(pieces) >= 1

    # points inside P \ Q
    @test any([0.25, 0.25] in R for R in pieces)
    @test any([1.75, 1.75] in R for R in pieces)

    # point inside Q should not belong to any piece
    @test !any([1.0, 1.0] in R for R in pieces)
end

@testset "SemiLinearSet difference decomposition" begin
    P1 = hbox([0.0, 0.0], [2.0, 2.0])
    P2 = hbox([3.0, 3.0], [4.0, 4.0])
    Q = hbox([0.5, 0.5], [1.5, 1.5])

    S = UT.semilinear_set([P1, P2])

    D1 = UT.semilinear_set(UT.set_difference_decompose(S, UT.semilinear_set(Q)))
    @test D1 isa UT.SemiLinearSet
    @test [0.25, 0.25] ∈ D1
    @test [3.5, 3.5] ∈ D1
    @test [1.0, 1.0] ∉ D1

    D2 = UT.semilinear_set(UT.set_difference_decompose(S, Q))
    @test D2 isa UT.SemiLinearSet
    @test [0.25, 0.25] ∈ D2
    @test [3.5, 3.5] ∈ D2
    @test [1.0, 1.0] ∉ D2
end

@testset "disjointify and semilinear volume" begin
    P1 = hbox([0.0, 0.0], [2.0, 1.0])   # area = 2
    P2 = hbox([1.0, 0.0], [3.0, 1.0])   # area = 2
    # overlap = [1,2]x[0,1], area = 1
    # union area = 3

    S = UT.semilinear_set([P1, P2])

    Sd = UT.disjointify(S)
    @test Sd isa UT.SemiLinearSet

    vd = UT.get_volume(Sd; backend = backend, assume_disjoint = true)
    v = UT.get_volume(S; backend = backend)

    @test isapprox(vd, 3.0; atol = 1e-8)
    @test isapprox(v, 3.0; atol = 1e-8)
end

@testset "preimage under linear map" begin
    P = hbox([0.0, 0.0], [1.0, 1.0])
    A = [2.0 0.0; 0.0 3.0]

    Ppre = UT.preimage_linear(P, A)
    @test Ppre isa Poly

    # Ax in [0,1]x[0,1]  => x in [0,1/2]x[0,1/3]
    @test [0.25, 0.2] ∈ Ppre
    @test [0.6, 0.2] ∉ Ppre
    @test [0.25, 0.4] ∉ Ppre

    S = UT.semilinear_set([P])
    Spre = UT.preimage_linear(S, A)
    @test Spre isa UT.SemiLinearSet
    @test length(Spre.array) == 1
    @test [0.25, 0.2] ∈ Spre
    @test [0.6, 0.2] ∉ Spre
end

@testset "preimage_linear_parts" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])
    S = UT.semilinear_set([P1, P2])

    A = Matrix{Float64}(LA.I, 2, 2)
    parts = UT.preimage_linear_parts(S, A)

    @test parts isa Vector
    @test length(parts) == 2
    @test [0.5, 0.5] ∈ parts[1] || [0.5, 0.5] ∈ parts[2]
    @test [2.5, 2.5] ∈ parts[1] || [2.5, 2.5] ∈ parts[2]
end

end
