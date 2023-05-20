module TestMain
    using Test
    Usz = 100
    Wsz = 0
    n_step = 2
    no_plot = true
    include("../../examples/example_SF_abst_2-simple.jl")
    
    @testset "state_trans" begin
        @test contr.data[1] == (18, 14)
        @test costBound ≈ 0.6757639538776234 rtol=1e-3
        @test costTrue <= costBound
    end
end
