module TestMain
    using Test
    Usz = 100
    Wsz = 0
    n_step = 2
    no_plot = true
    include("../../examples/example_SF_abst_2.jl")
    
    @testset "state_trans" begin
        @test contr.data[1] == (45,21)
        @test costBound ≈ 2.796543525044841 rtol=1e-3
        @test costTrue <= costTrue
    end
end
