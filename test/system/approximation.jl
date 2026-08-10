module TestKernelApproximations

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
using Random
import MathematicalSystems as MS

@testset "Simulation utilities" begin
    # Simple linear dynamics: ẋ = A*x + B*u
    A = @SMatrix [0.0 1.0; -1.0 0.0]
    B = @SMatrix [0.0; 1.0]
    dynamics(x, u) = A*x + B*u

    x0 = @SVector [1.0, 0.0]
    u0 = @SVector [0.2]
    tstep = 0.1

    x1 = ST.runge_kutta4(dynamics, x0, u0, tstep, 10)
    @test length(x1) == 2
    @test all(isfinite, x1)

    simmap = ST.simulate_control_map(dynamics; num_substeps = 7)
    x2 = simmap(x0, u0, tstep)
    @test x2 ≈ ST.runge_kutta4(dynamics, x0, u0, tstep, 7)

    discmap = ST.discretize_control_map(dynamics, tstep; num_substeps = 9)
    x3 = discmap(x0, u0)
    @test x3 ≈ ST.runge_kutta4(dynamics, x0, u0, tstep, 9)
end

@testset "Approximations: over/under maps" begin
    # ----------------------------
    # Make toy systems (discrete + continuous)
    # ----------------------------
    # discrete dynamics: x⁺ = x + u
    fd(x, u) = x + u

    # continuous dynamics: ẋ = u  (so x(t+h) = x + h*u)
    fc(x, u) = u

    n = 2
    m = 2

    X = UT.box([-10.0, -10.0], [10.0, 10.0])
    U = UT.box([-1.0, -1.0], [1.0, 1.0])

    sysD = MS.ConstrainedBlackBoxControlDiscreteSystem(fd, n, m, X, U)
    sysC = MS.ConstrainedBlackBoxControlContinuousSystem(fc, n, m, X, U)

    # Common test set and input
    rect = UT.box([-0.5, -0.5], [0.5, 0.5])
    u = @SVector [0.2, -0.1]
    tstep = 0.3

    # ----------------------------
    # DiscreteTimeOverApproximationMap
    # ----------------------------
    over_map_rect = (r, u) -> begin
        x = LazySets.center(r)
        Fx = fd(x, u)
        rad = LazySets.radius_hyperrectangle(r)
        UT.box(Fx - rad, Fx + rad)
    end
    Odisc = ST.DiscreteTimeOverApproximationMap(sysD, over_map_rect)
    @test ST.is_over_approximation(Odisc)
    @test ST.get_system(Odisc) === sysD

    outR = ST.get_over_approximation_map(Odisc)(rect, u)
    @test outR isa LazySets.AbstractHyperrectangle
    @test LazySets.center(outR) ≈ fd(LazySets.center(rect), u)

    # ----------------------------
    # ContinuousTimeOverApproximationMap + discretize
    # ----------------------------
    over_map_cont = (r, u, h) -> begin
        x = LazySets.center(r)
        Fx = x + h*u # exact for ẋ=u
        rad = LazySets.radius_hyperrectangle(r)
        UT.box(Fx - rad, Fx + rad)
    end
    Ocont = ST.ContinuousTimeOverApproximationMap(sysC, over_map_cont)
    @test ST.is_over_approximation(Ocont)

    outRc = ST.get_over_approximation_map(Ocont)(rect, u, tstep)
    @test LazySets.center(outRc) ≈ (LazySets.center(rect) + tstep*u)

    OcontD = ST.discretize(Ocont, tstep)
    outRc2 = ST.get_over_approximation_map(OcontD)(rect, u)
    @test LazySets.center(outRc2) ≈ LazySets.center(outRc)

    # ----------------------------
    # Growth bounds (discrete/continuous)
    # ----------------------------
    gb_disc = (r, u) -> abs.(r) .+ 0.1 .* abs.(u)  # simple monotone bound
    Gdisc = ST.DiscreteTimeGrowthBound(sysD, gb_disc)
    Rg = ST.get_over_approximation_map(Gdisc)(rect, u)
    @test Rg isa LazySets.AbstractHyperrectangle

    gb_cont = (r, u, h) -> abs.(r) .+ h .* 0.1 .* abs.(u)
    Gcont = ST.ContinuousTimeGrowthBound(sysC, gb_cont)
    Rgc = ST.get_over_approximation_map(Gcont)(rect, u, tstep)
    @test Rgc isa LazySets.AbstractHyperrectangle

    GcontD = ST.discretize(Gcont, tstep)
    Rgc2 = ST.get_over_approximation_map(GcontD)(rect, u)
    @test Rgc2 isa LazySets.AbstractHyperrectangle

    # ----------------------------
    # Linearized (discrete/continuous)
    # ----------------------------
    # linsys_map must return (Fx, DFx)
    # for x⁺ = x + u: DFx = I
    linsys_d = (x, dx, u) -> (fd(x, u), @SMatrix [1.0 0.0; 0.0 1.0])
    err_d = (e, u) -> 0.01 .* ones(SVector{2, Float64})
    Ldisc = ST.DiscreteTimeLinearized(sysD, linsys_d, err_d)
    Rl = ST.get_over_approximation_map(Ldisc)(rect, u)
    @test Rl isa LazySets.AbstractHyperrectangle
    @test LazySets.center(Rl) ≈ fd(LazySets.center(rect), u)

    linsys_c = (x, dx, u, h) -> (x + h*u, @SMatrix [1.0 0.0; 0.0 1.0])
    err_c = (e, u, h) -> 0.01 .* ones(SVector{2, Float64})
    Lcont = ST.ContinuousTimeLinearized(sysC, linsys_c, err_c)
    Rlc = ST.get_over_approximation_map(Lcont)(rect, u, tstep)
    @test LazySets.center(Rlc) ≈ (LazySets.center(rect) + tstep*u)

    LcontD = ST.discretize(Lcont, tstep)
    Rlc2 = ST.get_over_approximation_map(LcontD)(rect, u)
    @test LazySets.center(Rlc2) ≈ LazySets.center(Rlc)

    # ----------------------------
    # Under-approximations: centered
    # ----------------------------
    Cdisc = ST.DiscreteTimeCenteredSimulation(sysD)
    pts = ST.get_under_approximation_map(Cdisc)(rect, u)
    @test length(pts) == 1
    @test pts[1] ≈ fd(LazySets.center(rect), u)

    Ccont = ST.ContinuousTimeCenteredSimulation(sysC)
    pts2 = ST.get_under_approximation_map(Ccont)(rect, u, tstep)
    @test length(pts2) == 1
    @test pts2[1] ≈ (LazySets.center(rect) + tstep*u)

    CcontD = ST.discretize(Ccont, tstep)
    pts3 = ST.get_under_approximation_map(CcontD)(rect, u)
    @test pts3[1] ≈ pts2[1]

    # ----------------------------
    # Under-approximations: random
    # ----------------------------
    Random.seed!(123)
    Rdisc = ST.DiscreteTimeRandomSimulation(sysD, 20)
    ptsR = ST.get_under_approximation_map(Rdisc)(rect, u)
    @test length(ptsR) == 20
    @test all(p -> length(p) == 2, ptsR)

    Random.seed!(123)
    Rcont = ST.ContinuousTimeRandomSimulation(sysC, 15)
    ptsRc = ST.get_under_approximation_map(Rcont)(rect, u, tstep)
    @test length(ptsRc) == 15

    RcontD = ST.discretize(Rcont, tstep)
    ptsRd = ST.get_under_approximation_map(RcontD)(rect, u)
    @test length(ptsRd) == 15

    # Reproducibility with an explicit rng (also survives discretize)
    R1 = ST.DiscreteTimeRandomSimulation(sysD, 10, Random.MersenneTwister(42))
    R2 = ST.DiscreteTimeRandomSimulation(sysD, 10, Random.MersenneTwister(42))
    @test ST.get_under_approximation_map(R1)(rect, u) ==
          ST.get_under_approximation_map(R2)(rect, u)

    C1 = ST.ContinuousTimeRandomSimulation(sysC, 10, Random.MersenneTwister(7))
    C2 = ST.discretize(
        ST.ContinuousTimeRandomSimulation(sysC, 10, Random.MersenneTwister(7)),
        tstep,
    )
    @test ST.get_under_approximation_map(C1)(rect, u, tstep) ==
          ST.get_under_approximation_map(C2)(rect, u)
end

@testset "input_cache / reach_set hooks" begin
    fd(x, u) = x + u
    fc(x, u) = u
    X = UT.box([-10.0, -10.0], [10.0, 10.0])
    U = UT.box([-1.0, -1.0], [1.0, 1.0])
    sysD = MS.ConstrainedBlackBoxControlDiscreteSystem(fd, 2, 2, X, U)

    rect = UT.box([-0.5, -0.25], [0.5, 0.25])
    r = LazySets.radius_hyperrectangle(rect)
    u = @SVector [0.2, -0.1]

    # Growth bound: hook path ≡ generic map
    gb = (r, u) -> abs.(r) .+ 0.1 .* abs.(u)
    G = ST.DiscreteTimeGrowthBound(sysD, gb)
    cache = ST.input_cache(G, r, u)
    @test cache == gb(r, u)
    out_hook = ST.reach_set(G, rect, u, cache)
    out_map = ST.get_over_approximation_map(G)(rect, u)
    @test LazySets.center(out_hook) == LazySets.center(out_map)
    @test LazySets.radius_hyperrectangle(out_hook) ==
          LazySets.radius_hyperrectangle(out_map)

    # Linearized: hook path ≡ generic map
    linsys_d = (x, dx, u) -> (fd(x, u), @SMatrix [1.0 0.0; 0.0 1.0])
    err_d = (e, u) -> 0.01 .* ones(SVector{2, Float64})
    L = ST.DiscreteTimeLinearized(sysD, linsys_d, err_d)
    cacheL = ST.input_cache(L, r, u)
    out_hookL = ST.reach_set(L, rect, u, cacheL)
    out_mapL = ST.get_over_approximation_map(L)(rect, u)
    @test LazySets.center(out_hookL) == LazySets.center(out_mapL)
    @test LazySets.radius_hyperrectangle(out_hookL) ==
          LazySets.radius_hyperrectangle(out_mapL)

    # Default hooks: cache is nothing, reach_set falls back to the generic map
    over_map = (rc, uu) -> begin
        x = LazySets.center(rc)
        rad = LazySets.radius_hyperrectangle(rc)
        UT.box(fd(x, uu) - rad, fd(x, uu) + rad)
    end
    O = ST.DiscreteTimeOverApproximationMap(sysD, over_map)
    @test ST.input_cache(O, r, u) === nothing
    out_def = ST.reach_set(O, rect, u, nothing)
    @test LazySets.center(out_def) ≈ fd(LazySets.center(rect), u)
end

end # module
