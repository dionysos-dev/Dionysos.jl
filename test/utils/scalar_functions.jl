module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "ScalarFunctions" begin
    x_1 = 1.0
    f_0 = UT.ZeroFunction()
    f_1 = UT.ConstantFunction(1.0)
    f_2 = UT.ConstantFunction(2.0)

    @test UT.function_value(f_1, x_1) == 1.0
    @test UT.function_value(f_0 + f_1, x_1) == 1.0
    @test UT.function_value(f_1 + f_2, x_1) == 3.0

    x_2 = Vector{Float64}([1.0, 1.0])
    f_4 = UT.AffineFunction([1.0, 2.0], 3.0)
    @test UT.function_value(f_4, x_2) == 6.0

    f_4 = UT.QuadraticControlFunction(ones(Float64, (2, 2)))
    @test UT.function_value(f_4, x_2) == 4.0
end

println("End test")
end
