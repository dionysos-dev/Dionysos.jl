module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import Polyhedra

@testset "ScalarFunctions (callable contract)" begin
    x_1 = 1.0
    f_0 = UT.ZeroFunction()
    f_1 = UT.ConstantFunction(1.0)
    f_2 = UT.ConstantFunction(2.0)

    @test f_0(x_1) == 0.0
    @test f_0(x_1, [0.5]) == 0.0
    @test f_1(x_1) == 1.0
    @test (f_0 + f_1)(x_1) == 1.0
    @test (f_1 + f_2)(x_1) == 3.0

    x_2 = Vector{Float64}([1.0, 1.0])
    f_a = UT.AffineFunction([1.0, 2.0], 3.0)
    @test f_a(x_2) == 6.0
    @test f_a isa UT.ScalarFunction

    f_q = UT.QuadraticFunction(ones(Float64, (2, 2)))
    @test f_q(x_2) == 4.0

    u = [0.5]
    f_c = UT.ConstantControlFunction(7.0)
    @test f_c(x_2, u) == 7.0
    @test f_c isa UT.ScalarControlFunction

    f_qsc = UT.QuadraticStateControlFunction(
        ones(Float64, 2, 2),
        ones(Float64, 1, 1),
        zeros(2, 1),
        zeros(2),
        zeros(1),
        1.0,
    )
    # x'Qx + u'Ru + v = 4 + 0.25 + 1
    @test f_qsc(x_2, u) == 5.25
    @test UT.get_full_psd_matrix(f_qsc) == [
        1.0 1.0 0.0 0.0
        1.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]

    # any plain callable is a valid cost — the contract the solvers rely on
    user_cost(x, u) = sum(abs2, x) + sum(abs2, u)
    @test user_cost(x_2, u) == 2.25

    # black-box wrappers put user functions into the type hierarchy
    f_bb = UT.BlackBoxFunction(x -> sum(abs2, x))
    @test f_bb isa UT.ScalarFunction
    @test f_bb(x_2) == 2.0
    f_bbc = UT.BlackBoxControlFunction(user_cost)
    @test f_bbc isa UT.ScalarControlFunction
    @test f_bbc(x_2, u) == 2.25
end

@testset "AffineFunction approximate equality" begin
    f = UT.AffineFunction([1.0, 2.0], 3.0)
    @test isapprox(f, UT.AffineFunction([1.0, 2.0], 3.0))
    @test isapprox(f, UT.AffineFunction([1.0, 2.0 + 1e-12], 3.0))
    @test !isapprox(f, UT.AffineFunction([1.0, 2.5], 3.0))   # slope differs
    @test !isapprox(f, UT.AffineFunction([1.0, 2.0], 3.5))   # offset differs
end

@testset "PolyhedralFunction (max of pieces on a domain, ∞ outside)" begin
    # domain: the 1-D interval x ∈ [0, 2]
    domain = Polyhedra.HalfSpace([-1.0], 0.0) ∩ Polyhedra.HalfSpace([1.0], 2.0)
    pieces = [UT.AffineFunction([1.0], 0.0), UT.AffineFunction([-1.0], 2.0)]  # x and 2 − x
    p = UT.PolyhedralFunction(0.0, pieces, domain)

    @test p([0.5]) == 1.5   # max(0, 0.5, 1.5)
    @test p([1.0]) == 1.0   # max(0, 1.0, 1.0), the pieces meet
    @test p([3.0]) == Inf   # outside the domain

    # + ConstantFunction only lifts the lower bound (pieces and domain unchanged)
    shifted = UT.ConstantFunction(10.0) + p
    @test shifted isa UT.PolyhedralFunction
    @test shifted([0.5]) == 10.0  # the raised floor now dominates the pieces
    @test shifted([3.0]) == Inf   # still ∞ outside the domain

    # + ZeroFunction is the identity
    @test (UT.ZeroFunction() + p) === p

    # _inf is only defined for floats
    @test_throws ErrorException UT._inf(Int)
end

@testset "_cost_matrix (canonical [x; u; 1] PSD form)" begin
    # a raw matrix is assumed to already be in canonical form → identity
    M = [1.0 0.0; 0.0 2.0]
    @test UT._cost_matrix(M) === M

    # a QuadraticStateControlFunction converts through get_full_psd_matrix
    f = UT.QuadraticStateControlFunction(
        ones(2, 2),
        ones(1, 1),
        zeros(2, 1),
        zeros(2),
        zeros(1),
        1.0,
    )
    @test UT._cost_matrix(f) == UT.get_full_psd_matrix(f)
end

end
