module TestMain
using Test
no_plot = true
include("../../examples/example_SF_abst_2.jl")
@testset "state_trans" begin
    @test cost_bound ≈ 0.6250139513432214 rtol = 1e-3
    @test cost_true ≈ 0.36844089806471475 rtol = 1e-3
    @test cost_true <= cost_bound
end
end
