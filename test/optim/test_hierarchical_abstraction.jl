module TestMain
no_plot = true
include("../../examples/example_hierarchical_abstraction.jl")
@testset "hierarchical-abstraction" begin
    @test x_traj[end] ∈ concrete_problem.target_set
    @test cost === 38.0
end

end
