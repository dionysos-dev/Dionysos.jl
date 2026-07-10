module TestTrajectories

using Test
using StaticArrays
import MathematicalSystems as MS
import HybridSystems

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

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

@testset "DiscreteTrajectory" begin
    sys = DummyHybridSystem()

    traj0 = ST.DiscreteTrajectory{DummyTransition{Int}}(1)
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
