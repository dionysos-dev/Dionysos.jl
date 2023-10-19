module TestMain
include("../../utils/example_hierarchical_abstraction.jl")
@testset "hierarchical-abstraction" begin
    @test ST.get_state(cost_control_trajectory, ST.length(cost_control_trajectory)) ∈
          concrete_problem.target_set
    @test cost === 38.0
    @test isa(fig1, Plots.Plot{Plots.GRBackend})
    @test isa(fig2, Plots.Plot{Plots.GRBackend})
end

end
