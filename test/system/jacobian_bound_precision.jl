module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Symbolics
import MathematicalSystems as MS
import HybridSystems
import LazySets
import Random

# `compute_jacobian_bound` offers three levels, trading tightness against how often the bound
# is recomputed. All three are rigorous (interval arithmetic, not sampling); they differ only
# in how much of the state/input space each evaluation has to cover at once.
#
# The pendulum is the clean case: ∂f₂/∂x₁ = -(g/l)cos(x₁), so a bound over a full turn must be
# g/l, while a bound over a narrow slice of angle can be far smaller. That ordering is the
# property worth pinning — the absolute numbers are free to move.

const g, l = 9.81, 1.0
pendulum(X) = MS.ConstrainedBlackBoxControlContinuousSystem(
    (x, u) -> SVector(x[2], -(g / l) * sin(x[1]) + u[1]),
    2,
    1,
    X,
    UT.box(SVector(-3.0), SVector(3.0)),
)

@testset "compute_jacobian_bound precision levels" begin
    sys = pendulum(UT.box(SVector(-π, -5.0), SVector(π, 5.0)))

    Lglobal = ST.compute_jacobian_bound(sys; precision = ST.GLOBAL_BOUND)
    Linput = ST.compute_jacobian_bound(sys; precision = ST.INPUT_BOUND)
    Lregion = ST.compute_jacobian_bound(sys; precision = ST.REGIONWISE_BOUND, nsplit = 8)

    @test Lregion isa ST.RegionwiseBound
    @test !(Linput isa ST.RegionwiseBound)
    @test Lregion.nregions == 8^2

    u = SVector(0.0)

    # Over the whole turn the cosine reaches 1, so every level must give g/l there.
    @test Matrix(Lglobal(u))[2, 1] ≈ g / l rtol = 1e-6
    @test Matrix(Linput(u))[2, 1] ≈ g / l rtol = 1e-6

    at(x) = Matrix(Lregion.bound(Lregion.region_of(x), u))
    # Near x₁ = π/2 the cosine is small, and only the regionwise level can see that.
    @test at(SVector(π / 2, 0.0))[2, 1] < g / l
    # ... while at the bottom, where |cos| is largest, it must not undercut the truth.
    @test at(SVector(0.0, 0.0))[2, 1] >= abs((g / l) * cos(0.0)) - 1e-6
end

@testset "regionwise bound stays sound" begin
    # Whatever the level, the bound must dominate the real Jacobian at every sampled point of
    # the region it is claimed for — that is the whole contract.
    sys = pendulum(UT.box(SVector(-π, -5.0), SVector(π, 5.0)))
    L = ST.compute_jacobian_bound(sys; precision = ST.REGIONWISE_BOUND, nsplit = 8)
    rng = Random.MersenneTwister(1)
    worst = -Inf
    for _ in 1:2000
        x = SVector(-π + 2π * rand(rng), -5.0 + 10.0 * rand(rng))
        u = SVector(-3.0 + 6.0 * rand(rng))
        M = Matrix(L.bound(L.region_of(x), u))
        worst = max(worst, abs(-(g / l) * cos(x[1])) - M[2, 1])
        worst = max(worst, 1.0 - M[1, 2])
    end
    @test worst <= 1e-6
end

@testset "a regionwise bound reaches the abstraction, and stays hoisted" begin
    # The efficiency property: the radius must be integrated once per *region*, not once per
    # cell. `input_cache` therefore still returns a table — one entry per region — and the
    # per-cell work is a lookup.
    sys = pendulum(UT.box(SVector(-π, -5.0), SVector(π, 5.0)))
    L = ST.compute_jacobian_bound(sys; precision = ST.REGIONWISE_BOUND, nsplit = 8)
    disc = ST.discretize(ST.ContinuousTimeGrowthBound(sys; jacobian_bound = L), 0.1)

    cache = ST.input_cache(disc, SVector(0.05, 0.05), SVector(0.0))
    @test cache isa AbstractVector
    @test length(cache) == 8^2

    # Two cells of the same size at different angles must get different radii.
    cell(c) = LazySets.Hyperrectangle(c, SVector(0.05, 0.05))
    over = ST.get_over_approximation_map(disc)
    r_bottom = LazySets.radius_hyperrectangle(over(cell(SVector(0.0, 0.0)), SVector(0.0)))
    r_quarter =
        LazySets.radius_hyperrectangle(over(cell(SVector(π / 2, 0.0)), SVector(0.0)))
    @test r_bottom != r_quarter
    # Tighter where the dynamics are milder — the whole point of splitting.
    @test r_quarter[2] < r_bottom[2]

    # An input-only bound still takes the plain hoisted path.
    Li = ST.compute_jacobian_bound(sys; precision = ST.INPUT_BOUND)
    disc_i = ST.discretize(ST.ContinuousTimeGrowthBound(sys; jacobian_bound = Li), 0.1)
    @test ST.input_cache(disc_i, SVector(0.05, 0.05), SVector(0.0)) isa AbstractVector
end

@testset "precision is selectable through MOI" begin
    # The knob is useless if it cannot be reached from the interface every solver is driven
    # through, so drive it that way: set the attribute, solve, and check the abstraction that
    # comes out differs from the default one.
    import MathOptInterface as MOI

    function abstraction(; precision = nothing, nsplit = nothing)
        X = UT.box(SVector(-π, -5.0), SVector(π, 5.0))
        problem = PR.AlternatingSimulationProblem(pendulum(X), X)
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("state_grid"),
            MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3)),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(0.0), SVector(1.0)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), 0.1)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.GROWTH,
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
        precision === nothing ||
            MOI.set(opt, MOI.RawOptimizerAttribute("jacobian_bound_precision"), precision)
        nsplit === nothing ||
            MOI.set(opt, MOI.RawOptimizerAttribute("jacobian_bound_nsplit"), nsplit)
        MOI.optimize!(opt)
        return MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
    end

    ntrans(sys) = HybridSystems.ntransitions(SY.get_automaton(sys))

    # A regionwise bound is tighter, so it must not produce *more* spurious transitions than
    # the whole-domain one — that is the abstraction quality the precision buys.
    coarse = ntrans(abstraction(; precision = ST.GLOBAL_BOUND))
    fine = ntrans(abstraction(; precision = ST.REGIONWISE_BOUND, nsplit = 8))
    @test fine <= coarse
    # And the default is still reachable and sane.
    @test ntrans(abstraction()) > 0
end

end # module TestMain
