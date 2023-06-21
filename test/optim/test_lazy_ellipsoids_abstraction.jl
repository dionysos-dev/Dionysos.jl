module TestMain
no_plot = true
include("../../examples/SFRRT/example_lazy_ellipsoid_1.jl")
@testset "lazy-ellipsoids-abstraction" begin
    @test cost_true <= cost_bound
end

end
