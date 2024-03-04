using Test
using Dionysos
newton_method = Dionysos.Utils.newton_method

# Define the test function, its derivative, and second derivative
@testset "EllipsoidBasics" begin
    f(x) = x^2 - 4
    df(x) = 2x
    ddf(x) = 2

    # Test parameters
    interval = [1, 3]
    x0 = 1.5
    ϵ = 1e-6

    # Call the newton_method function with the test parameters
    result = newton_method(
        f,
        df,
        ddf;
        interval = interval,
        x0 = x0,
        ϵ = ϵ,
        verbose = false,
        stopIfNegative = false,
    )
    # Check the result
    @test isapprox(result[1], -3.0; atol = ϵ) && interval[1] ≤ result[2] ≤ interval[2]

    result = newton_method(
        f,
        df,
        ddf;
        interval = [0, 4],
        x0 = x0,
        ϵ = ϵ,
        verbose = false,
        stopIfNegative = false,
    )
    @test isapprox(result[1], -4.0; atol = ϵ) && 0 ≤ result[2] ≤ 4
end
