module TestTrajectories

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS
import HybridSystems
import Suppressor
using Plots

# ------------------------------------------------------------
# Minimal controller types for tests
# ------------------------------------------------------------
struct DummyMap{F}
    h::F
end

# ------------------------------------------------------------
# Minimal "hybrid system" for last_mode tests
# We only need HybridSystems.target(system, transition)
# ------------------------------------------------------------
struct DummyHybridSystem end
struct DummyTransition{Q}
    tgt::Q
end
HybridSystems.target(::DummyHybridSystem, t::DummyTransition) = t.tgt

@testset "AutomatonPath" begin
    sys = DummyHybridSystem()

    traj0 = ST.AutomatonPath{DummyTransition{Int}}(1)
    @test traj0.q_0 == 1
    @test length(traj0) == 0
    @test ST.last_mode(sys, traj0) == 1

    t1 = DummyTransition(2)
    traj1 = ST.append(traj0, t1)
    @test length(traj1) == 1
    @test traj1.transitions[1] == t1
    @test ST.last_mode(sys, traj1) == 2

    t2 = DummyTransition(5)
    traj2 = ST.append(traj1, t2)
    @test length(traj2) == 2
    @test ST.last_mode(sys, traj2) == 5
end

@testset "Trajectory container + cost" begin
    xs = [@SVector([0.0, 0.0]), @SVector([1.0, 0.0]), @SVector([1.0, 2.0])]
    us = [@SVector([1.0, 0.0]), @SVector([0.0, 2.0])]
    traj = ST.Trajectory(xs; inputs = us)

    @test length(traj) == 3
    @test ST.states(traj) === xs
    @test ST.inputs(traj) === us
    @test ST.times(traj) === nothing
    @test ST.modes(traj) === nothing
    @test ST.memory(traj) === nothing

    c(x, u) = sum(abs, x) + sum(abs, u)
    out = ST.get_cost_trajectory(traj, c)
    @test out.cost isa Vector
    @test length(out.cost) == 2

    # manual total:
    c1 = c(xs[1], us[1])
    c2 = c(xs[2], us[2])
    @test out.total_cost ≈ (c1 + c2)

    # mismatched channel lengths are rejected at construction (3 states need 2 inputs)
    @test_throws DimensionMismatch ST.Trajectory(xs; inputs = [@SVector([0.0, 0.0])])
    @test_throws DimensionMismatch ST.Trajectory(xs; modes = [1, 1])

    # a state-only trajectory has no input channel
    stateonly = ST.Trajectory(xs)
    @test ST.inputs(stateonly) === nothing
    @test_throws ErrorException ST.get_cost_trajectory(stateonly, c)
end

@testset "Closed-loop engine" begin
    # discrete dynamics x⁺ = x + u (broadcast tolerates scalar or vector u)
    f(x, u) = x .+ u
    sys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 1, 1, nothing, nothing)

    # static controller under test: u = -x/2
    ctrl = ST.AffineController(MS.AffineMap(hcat(-0.5), [0.0]))

    traj = ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 5)
    @test traj isa ST.Trajectory
    @test ST.memory(traj) === nothing
    @test length(ST.states(traj)) == 6
    @test length(ST.inputs(traj)) == 5
    @test ST.states(traj)[2] ≈ [0.5]                # x + (-0.5x)
    @test eltype(ST.inputs(traj)) == Vector{Float64}  # inputs typed from first control

    # channels of a plain continuous rollout
    traj3 = ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 3)
    @test length(ST.states(traj3)) == 4 && length(ST.inputs(traj3)) == 3
    @test ST.modes(traj3) === nothing && ST.times(traj3) === nothing

    # stopping condition
    traj_stop =
        ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 50; stopping = x -> x[1] < 0.1)
    @test ST.states(traj_stop)[end][1] < 0.1

    # dynamic controller: memory counts steps and cuts off after 2 controls
    dyn = ST.DiscreteDynamicController(
        0,
        ST.PredicateDomain(memx -> memx[1] < 2),
        (mem, y) -> mem + 1,
        (mem, y) -> mem < 2 ? [0.5] : nothing,
        false,
    )
    trajd = ST.get_closed_loop_trajectory(sys, dyn, [0.0], 10)
    @test length(ST.inputs(trajd)) == 2             # third output_control is undefined
    @test ST.memory(trajd) == [0, 1, 2]
end

@testset "Closed-loop engine: early-stop and verbose branches" begin
    # x⁺ = x + u; a static feedback that would otherwise run forever.
    f(x, u) = x .+ u
    sys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 1, 1, nothing, nothing)
    ctrl = ST.AffineController(MS.AffineMap(hcat(-0.5), [0.0]))

    # Each scenario trips a distinct early-stop branch. `verbose = true` exercises the
    # log-emitting paths; `@suppress` keeps their output out of the test log.

    # (1) `trajectory_success` already holds for the initial one-state trajectory.
    t = Suppressor.@suppress ST.get_closed_loop_trajectory(
        sys,
        ctrl,
        [1.0],
        5;
        trajectory_success = _ -> true,
        verbose = true,
    )
    @test length(ST.states(t)) == 1

    # (2) the controller yields a valid input but the system has no successor.
    t = Suppressor.@suppress ST.get_closed_loop_trajectory(
        sys,
        ctrl,
        [1.0],
        5;
        f_map_override = (x, u) -> nothing,
        verbose = true,
    )
    @test length(ST.states(t)) == 1

    # (3) the successor state is non-finite.
    t = Suppressor.@suppress ST.get_closed_loop_trajectory(
        sys,
        ctrl,
        [1.0],
        5;
        f_map_override = (x, u) -> SVector(NaN),
        verbose = true,
    )
    @test length(ST.states(t)) == 1

    # (4) a dynamic controller whose memory update aborts the run: the output map is
    #     always admissible, but the state map returns `nothing` on the first step.
    dyn = ST.DiscreteDynamicController(
        0,
        ST.PredicateDomain(_ -> true),
        (mem, y) -> nothing,
        (mem, y) -> [0.5],
        false,
    )
    t = Suppressor.@suppress ST.get_closed_loop_trajectory(
        sys,
        dyn,
        [0.0],
        5;
        verbose = true,
    )
    @test length(ST.states(t)) == 1
    @test isempty(ST.inputs(t))
end

@testset "Trajectory channels + accessors" begin
    # explicit timed/hybrid channels round-trip through the accessors
    xs = [@SVector([0.0]), @SVector([1.0]), @SVector([2.0])]
    us = [@SVector([1.0]), @SVector([1.0])]
    ts = [0.0, 0.3, 0.6]
    ks = [1, 1, 2]
    htraj = ST.Trajectory(xs; inputs = us, times = ts, modes = ks)
    @test ST.states(htraj) === xs
    @test ST.inputs(htraj) === us
    @test ST.times(htraj) == ts
    @test ST.modes(htraj) == ks
    @test length(htraj) == 3
end

@testset "wrap_coord + wrapper" begin
    x = @SVector [3.2, -1.0, 7.5]
    periodic_dims = SVector(1, 3)
    periods = SVector(4.0, 10.0)
    start = SVector(0.0, 5.0)   # dim1 in [0,4), dim3 in [5,15)

    w = UT.wrap_coord(x, periodic_dims, periods; start = start)

    # dim2 unchanged
    @test w[2] == x[2]

    # wrapped ranges
    @test 0.0 <= w[1] < 4.0
    @test 5.0 <= w[3] < 15.0

    # wrapper function
    wrapfun = UT.get_periodic_wrapper(periodic_dims, periods; start = start)
    @test wrapfun(x) == w
end

# Asserted through the full Plots pipeline (`plot` → `plt.series_list`); see
# `test/utils/plotting.jl` for why that is the level that actually exercises a recipe.
@testset "Trajectory recipe" begin
    xs = [SVector(0.0, 0.0, 9.0), SVector(1.0, 0.5, 8.0), SVector(2.0, 1.0, 7.0)]
    traj = ST.Trajectory(xs)

    # One marker per state, plus one arrow per step.
    plt = plot(traj)
    @test length(plt.series_list) == 3 + 2
    @test count(s -> s[:seriestype] === :scatter, plt.series_list) == 3
    @test count(s -> s[:seriestype] === :path, plt.series_list) == 2

    # `with_arrows = false` must actually drop the arrows. It is spelled `with_arrows` and not
    # `arrows` because Plots treats `arrows` as an alias of its own `arrow` attribute and rewrites
    # it before the recipe runs — under the old name the option was silently ignored, and every
    # example page asking for an arrow-free trajectory still got arrows.
    bare = plot(traj; with_arrows = false)
    @test length(bare.series_list) == 3
    @test all(s -> s[:seriestype] === :scatter, bare.series_list)

    # Only the first point carries the label, so a trajectory is one legend entry, not 2n-1.
    @test [s[:label] for s in plot(traj; label = "closed loop").series_list] == ["closed loop"; fill("", 4)]
    @test [s[:label] for s in plot(traj).series_list] == fill("", 5)

    # `dims` chooses the plane, markers and arrows alike.
    proj = plot(traj; dims = [1, 3])
    @test proj.series_list[1][:y] == [9.0]
    @test proj.series_list[3][:y] == [9.0, 8.0]   # the first arrow

    # A single-state trajectory has no steps at all: one marker, no arrows.
    @test length(plot(ST.Trajectory([SVector(0.0, 0.0)])).series_list) == 1
end

end # module
