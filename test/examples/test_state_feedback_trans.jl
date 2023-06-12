module TestMain
    using Test
    Usz = 100
    Wsz = 0
    dt = 0.01
    n_step = 2
    simple = true
    no_plot = true
    include("../../examples/example_SF_abst_2.jl")
    @testset "state_trans" begin
        @test cost_bound ≈ 0.7384349757233135 rtol=1e-3
        @test cost_true ≈ 0.5691152071561884 rtol=1e-3
        @test cost_true <= cost_bound
    end
end
