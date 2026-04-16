
module TestMain

using Test
using Dionysos
using LazySets
using Polyhedra
using CDDLib
import LinearAlgebra as LA

const DI = Dionysos
const UT = DI.Utils

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

    C = UT.clean_poly(P)
    @test C isa Poly

    @test UT.dim(P) == 2
    @test UT.is_nonempty_set(P) == true

    E = LazySets.EmptySet{Float64}(2)
    @test UT.is_nonempty_set(E) == false
end

@testset "Poly intersection" begin
    P = hbox([0.0, 0.0], [1.0, 1.0])
    Q = hbox([0.5, 0.5], [2.0, 2.0])

    I = UT.set_intersection(P, Q)
    @test I isa Poly
    @test UT.is_nonempty_set(I)

    @test [0.75, 0.75] ∈ I
    @test [0.25, 0.25] ∉ I
    @test [1.5, 1.5] ∉ I
end

@testset "Poly volume" begin
    P = hbox([0.0, 0.0], [1.0, 2.0])
    v = UT.get_volume(P; backend = backend)
    @test isapprox(v, 2.0; atol = 1e-8)

    @test_throws Exception UT.get_volume(P)
end

@testset "SemiLinearSet constructors and basics" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])

    S0 = UT.SemiLinearSet()
    @test length(S0) == 0
    @test isempty(S0)

    S1 = UT.SemiLinearSet(P1)
    @test length(S1) == 1
    @test !isempty(S1)
    @test UT.dim(S1) == 2

    S2 = UT.SemiLinearSet([P1, P2])
    @test length(S2) == 2
    @test !isempty(S2)

    @test [0.5, 0.5] ∈ S2
    @test [2.5, 2.5] ∈ S2
    @test [1.5, 1.5] ∉ S2

    io = IOBuffer()
    show(io, S2)
    shown = String(take!(io))
    @test occursin("SemiLinearSet(2 parts)", shown)
end

@testset "normalize_semilinear" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])

    # Empty polytope in H-representation: x <= 0 and x >= 1 simultaneously
    E = Poly([HS([1.0, 0.0], 0.0), HS([-1.0, 0.0], -1.0)])

    S = UT.SemiLinearSet([P1, E])
    Sn = UT.normalize_semilinear(S)

    @test length(Sn) == 1
    @test [0.5, 0.5] ∈ Sn
    @test [2.0, 2.0] ∉ Sn
end

@testset "SemiLinearSet nonemptiness" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])

    E = Poly([HS([1.0, 0.0], 0.0), HS([-1.0, 0.0], -1.0)])

    S1 = UT.SemiLinearSet([E])
    S2 = UT.SemiLinearSet([E, P1])

    @test UT.is_nonempty_set(S1) == false
    @test UT.is_nonempty_set(S2) == true
end

@testset "SemiLinearSet intersection" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])
    Q = hbox([0.5, 0.5], [2.5, 2.5])

    S = UT.SemiLinearSet([P1, P2])

    I1 = UT.set_intersection(S, UT.SemiLinearSet(Q))
    @test I1 isa UT.SemiLinearSet
    @test length(I1) == 2
    @test [0.75, 0.75] ∈ I1
    @test [2.25, 2.25] ∈ I1
    @test [1.5, 1.5] ∉ I1

    I2 = UT.set_intersection(S, Q)
    @test I2 isa UT.SemiLinearSet
    @test [0.75, 0.75] ∈ I2
    @test [2.25, 2.25] ∈ I2
    @test [1.5, 1.5] ∉ I2

    I3 = UT.set_intersection(Q, S)
    @test I3 isa UT.SemiLinearSet
    @test [0.75, 0.75] ∈ I3
    @test [2.25, 2.25] ∈ I3
end

@testset "Poly difference decomposition" begin
    P = hbox([0.0, 0.0], [2.0, 2.0])
    Q = hbox([0.5, 0.5], [1.5, 1.5])

    pieces = UT.set_difference_decompose(P, Q)
    @test pieces isa Vector{Poly}
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

    S = UT.SemiLinearSet([P1, P2])

    D1 = UT.set_difference_decompose(S, UT.SemiLinearSet(Q))
    @test D1 isa UT.SemiLinearSet
    @test [0.25, 0.25] ∈ D1
    @test [3.5, 3.5] ∈ D1
    @test [1.0, 1.0] ∉ D1

    D2 = UT.set_difference_decompose(S, Q)
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

    S = UT.SemiLinearSet([P1, P2])

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

    S = UT.SemiLinearSet([P])
    Spre = UT.preimage_linear(S, A)
    @test Spre isa UT.SemiLinearSet
    @test length(Spre) == 1
    @test [0.25, 0.2] ∈ Spre
    @test [0.6, 0.2] ∉ Spre
end

@testset "preimage_linear_parts" begin
    P1 = hbox([0.0, 0.0], [1.0, 1.0])
    P2 = hbox([2.0, 2.0], [3.0, 3.0])
    S = UT.SemiLinearSet([P1, P2])

    A = Matrix{Float64}(LA.I, 2, 2)
    parts = UT.preimage_linear_parts(S, A)

    @test parts isa Vector
    @test length(parts) == 2
    @test [0.5, 0.5] ∈ parts[1] || [0.5, 0.5] ∈ parts[2]
    @test [2.5, 2.5] ∈ parts[1] || [2.5, 2.5] ∈ parts[2]
end

println("End test")
end
