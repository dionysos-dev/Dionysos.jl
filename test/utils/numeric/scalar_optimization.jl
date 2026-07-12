module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "golden_section_search" begin
    test_function(x) = x^2 - 4
    ϵ = 1e-6

    result = UT.golden_section_search(test_function; interval = [0, 4], δ = ϵ)
    @test isapprox(result[1], -4.0; atol = ϵ) && 0 ≤ result[2] ≤ 4

    result = UT.golden_section_search(test_function; interval = [1, 3], δ = ϵ)
    @test isapprox(result[1], -3.0; atol = ϵ) && 1 ≤ result[2] ≤ 3

    # stopIfNegative bails out before reaching the minimum, so the returned value
    # is negative but not the true minimum.
    result = UT.golden_section_search(
        test_function;
        interval = [-5, 4],
        δ = ϵ,
        stopIfNegative = true,
    )
    @test (-5 ≤ result[2] ≤ 4) && abs(result[1] + 4.0) > ϵ && result[1] < 0
end

@testset "newton_method" begin
    f(x) = x^2 - 4
    df(x) = 2x
    ddf(x) = 2
    ϵ = 1e-6

    result = UT.newton_method(f, df, ddf; interval = [1, 3], x0 = 1.5, ϵ = ϵ)
    @test isapprox(result[1], -3.0; atol = ϵ) && 1 ≤ result[2] ≤ 3

    result = UT.newton_method(f, df, ddf; interval = [0, 4], x0 = 1.5, ϵ = ϵ)
    @test isapprox(result[1], -4.0; atol = ϵ) && 0 ≤ result[2] ≤ 4
end

@testset "derivative_bisection" begin
    # convex f with a minimum at x = 1, f(1) = 2
    f(x) = (x - 1)^2 + 2
    df(x) = 2 * (x - 1)
    ddf(x) = 2
    δ = 1e-8

    result = UT.derivative_bisection(f, df, ddf; interval = [-5, 5], δ = δ)
    @test isapprox(result[2], 1.0; atol = 1e-4)
    @test isapprox(result[1], 2.0; atol = 1e-4)
end

end  # module TestMain
