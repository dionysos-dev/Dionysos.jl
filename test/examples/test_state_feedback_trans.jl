module TestMain
    using Test

    include("../../examples/example_SF_abst_2.jl")

    @testset "state_trans" begin
        @test contr.data[1] == (18,344)
        @test costBound ≈ 2.052186644292235 rtol=1e-3
    end
end