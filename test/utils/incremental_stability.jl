module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "Grid sizing from a Lyapunov matrix and per-step contraction" begin
    P = [2.0 0.0; 0.0 8.0]
    γ = 0.9          # per-step contraction of V(x, y) = ‖x - y‖_P
    ε = 0.1

    η = UT.grid_step_from_lyapunov(P, ε, γ)
    @test η > 0

    # η ↦ ε is the exact inverse of ε ↦ η.
    @test isapprox(UT.precision_from_grid_step(P, η, γ), ε; rtol = 1e-10)

    # Coarser precision ⇒ coarser grid.
    @test UT.grid_step_from_lyapunov(P, 2 * ε, γ) > η
    # Stronger contraction (smaller γ) ⇒ more slack per step ⇒ coarser grid allowed.
    @test UT.grid_step_from_lyapunov(P, ε, 0.5) > η
    # Worse-conditioned metric ⇒ finer grid required.
    P_ill = [2.0 0.0; 0.0 200.0]
    @test UT.grid_step_from_lyapunov(P_ill, ε, γ) < η
end

@testset "Non-contractive γ is rejected" begin
    P = [1.0 0.0; 0.0 1.0]
    @test_throws ErrorException UT.grid_step_from_lyapunov(P, 0.1, 1.0)
    @test_throws ErrorException UT.grid_step_from_lyapunov(P, 0.1, 1.5)
    @test_throws ErrorException UT.precision_from_grid_step(P, 0.1, 1.0)
end

end # module TestMain
