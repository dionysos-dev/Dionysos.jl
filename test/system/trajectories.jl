module TestTrajectories

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS
import HybridSystems

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
    xs = ST.Trajectory([@SVector([0.0, 0.0]), @SVector([1.0, 0.0]), @SVector([1.0, 2.0])])
    us = ST.Trajectory([@SVector([1.0, 0.0]), @SVector([0.0, 2.0])])

    @test length(xs) == 3
    @test length(us) == 2
    @test ST.get_elem(xs, 2) == @SVector([1.0, 0.0])
    @test ST.enum_elems(us) == us.seq

    c(x, u) = sum(abs, x) + sum(abs, u)
    out = ST.get_cost_trajectory(xs, us, c)
    @test out.c isa ST.Trajectory
    @test length(out.c) == 2

    # manual total:
    c1 = c(xs.seq[1], us.seq[1])
    c2 = c(xs.seq[2], us.seq[2])
    @test out.total_cost ≈ (c1 + c2)

    # mismatched lengths should assert
    bad_us = ST.Trajectory([@SVector([0.0, 0.0])])
    @test_throws AssertionError ST.get_cost_trajectory(xs, bad_us, c)
end

@testset "Closed-loop engine" begin
    # discrete dynamics x⁺ = x + u (broadcast tolerates scalar or vector u)
    f(x, u) = x .+ u
    sys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 1, 1, nothing, nothing)

    # static controller under test: u = -x/2
    ctrl = ST.AffineController(MS.AffineMap(hcat(-0.5), [0.0]))

    traj = ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 5)
    @test traj isa ST.ClosedLoopTrajectory
    @test traj.q === nothing
    @test length(traj.x) == 6
    @test length(traj.u) == 5
    @test traj.x.seq[2] ≈ [0.5]              # x + (-0.5x)
    @test eltype(traj.u.seq) == Vector{Float64}  # inputs typed from first control

    # destructuring (static: two elements)
    x_traj, u_traj = ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 3)
    @test length(x_traj) == 4 && length(u_traj) == 3

    # stopping condition
    traj_stop =
        ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 50; stopping = x -> x[1] < 0.1)
    @test traj_stop.x.seq[end][1] < 0.1

    # dynamic controller: memory counts steps and cuts off after 2 controls
    dyn = ST.DiscreteDynamicController(
        0,
        ST.PredicateDomain(memx -> memx[1] < 2),
        (mem, y) -> mem + 1,
        (mem, y) -> mem < 2 ? [0.5] : nothing,
        false,
    )
    trajd = ST.get_closed_loop_trajectory(sys, dyn, [0.0], 10)
    @test trajd.q isa ST.Trajectory
    @test length(trajd.u) == 2                # third output_control is undefined
    @test trajd.q.seq == [0, 1, 2]

    # dynamic destructuring: three elements
    xd, ud, qd = ST.get_closed_loop_trajectory(sys, dyn, [0.0], 10)
    @test qd.seq == trajd.q.seq
end

@testset "ClosedLoopTrajectory channels + accessors" begin
    f(x, u) = x .+ u
    sys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 1, 1, nothing, nothing)
    ctrl = ST.AffineController(MS.AffineMap(hcat(-0.5), [0.0]))

    # static rollout: state/input channels present, memory/time/mode absent
    traj = ST.get_closed_loop_trajectory(sys, ctrl, [1.0], 5)
    @test ST.states(traj) === traj.x.seq
    @test ST.inputs(traj) === traj.u.seq
    @test ST.memory(traj) === nothing
    @test ST.times(traj) === nothing
    @test ST.modes(traj) === nothing

    # explicit timed/hybrid channels round-trip through the accessors
    xs = ST.Trajectory([@SVector([0.0]), @SVector([1.0]), @SVector([2.0])])
    us = ST.Trajectory([@SVector([1.0]), @SVector([1.0])])
    ts = [0.0, 0.3, 0.6]
    ks = [1, 1, 2]
    htraj = ST.ClosedLoopTrajectory(xs, us; times = ts, modes = ks)
    @test ST.states(htraj) == xs.seq
    @test ST.inputs(htraj) == us.seq
    @test ST.times(htraj) == ts
    @test ST.modes(htraj) == ks

    # channels don't disturb the legacy destructuring contract
    x_traj, u_traj = htraj
    @test x_traj === xs && u_traj === us

    # dynamic controller memory still surfaces through `memory`
    dyn = ST.DiscreteDynamicController(
        0,
        ST.PredicateDomain(memx -> memx[1] < 2),
        (mem, y) -> mem + 1,
        (mem, y) -> mem < 2 ? [0.5] : nothing,
        false,
    )
    trajd = ST.get_closed_loop_trajectory(sys, dyn, [0.0], 10)
    @test ST.memory(trajd) == trajd.q.seq
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

end # module
