using Test
using Dionysos
bisection = Dionysos.Utils.bisection

# Define the test function, its derivative, and second derivative
@testset "Bissection" begin
    # Define a test function
    test_function(x) = x^2 - 4
    ϵ = 1e-6

    # Call the bisection function
    result = bisection(test_function; interval = [0, 4], δ = ϵ, verbose = false)
    @test isapprox(result[1], -4.0; atol = ϵ) && 0 ≤ result[2] ≤ 4

    result = bisection(test_function; interval = [0, 4], δ = ϵ, verbose = false)
    @test isapprox(result[1], -4.0; atol = ϵ) && 0 ≤ result[2] ≤ 4

    result = bisection(test_function; interval = [1, 3], δ = ϵ, verbose = false)
    @test isapprox(result[1], -3.0; atol = ϵ) && 0 ≤ result[2] ≤ 4

    result = bisection(
        test_function;
        interval = [-5, 4],
        δ = ϵ,
        verbose = false,
        stopIfNegative = true,
    )
    @test (-5 ≤ result[2] ≤ 4) && abs(result[1] + 4.0) > ϵ && result[1] < 0
end
