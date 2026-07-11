module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

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

println("End test")
end
