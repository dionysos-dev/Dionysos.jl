module TestMain

include("units_tests_LazyEllispoidAbstraction/example_lazy_ellipsoid.jl")
@testset "lazy-ellipsoids-abstraction" begin
    @test cost_true <= cost_bound
    @test isa(fig1, Plots.Plot{Plots.GRBackend})
    @test isa(fig2, Plots.Plot{Plots.GRBackend})
    @test isa(fig3, Plots.Plot{Plots.GRBackend})
end

end
