module TestMain
using Test
include("../../docs/src/examples/solvers/Hierarchical-abstraction.jl")
@testset "hierarchical-abstraction" begin
    @test x_traj.seq[end] ∈ concrete_problem.target_set
    @test cost === 38.0
    @test isa(fig1, Plots.Plot{Plots.GRBackend})
    @test isa(fig2, Plots.Plot{Plots.GRBackend})
end

end
