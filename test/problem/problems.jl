module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots
import MathematicalSystems
import LazySets

# Minimal dummy system: the problem plot recipes require `MS.stateset(system)`.
struct DummySystem{X}
    X::X
end
MathematicalSystems.stateset(s::DummySystem) = s.X

@testset "Problems" begin
    @testset "Infinity" begin
        @test Base.isfinite(4) == true
        @test Base.isfinite(-2.5) == true
        @test Base.isfinite(Inf) == false

        @test Base.isfinite(PR.Infinity()) == false
    end

    @testset "AlternatingSimulationProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
        sys = DummySystem(X)
        state_set =
            LazySets.Hyperrectangle(; low = SVector(-0.5, -0.5), high = SVector(0.5, 0.5))

        p = PR.AlternatingSimulationProblem(sys, state_set)

        @test p isa PR.ProblemType
        @test p isa PR.AbstractionProblem
        @test p.system === sys
        @test p.state_set == state_set

        # problems are immutable value objects
        state_set2 =
            LazySets.Hyperrectangle(; low = SVector(-0.2, -0.2), high = SVector(0.2, 0.2))
        @test_throws ErrorException p.state_set = state_set2
    end

    @testset "OptimalControlProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        sys = DummySystem(X)

        XI =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(-0.5, -0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(1.0, 1.0))

        state_cost = x -> 0.0
        transition_cost = (x, u) -> 1.0
        T = 5

        p = PR.OptimalControlProblem(sys, XI, XT, state_cost, transition_cost, T)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.target_set == XT
        @test p.state_cost === state_cost
        @test p.transition_cost === transition_cost
        @test p.time == T
        @test p.safe_set === nothing        # optional: no avoid part unless asked for

        XS = LazySets.Hyperrectangle(; low = SVector(-1.5, -1.5), high = SVector(1.5, 1.5))
        avoid = PR.OptimalControlProblem(
            sys,
            XI,
            XT,
            state_cost,
            transition_cost,
            T;
            safe_set = XS,
        )
        @test avoid.safe_set == XS
        # The five-argument form keeps defaulting the horizon, and still takes a safe set.
        @test PR.OptimalControlProblem(sys, XI, XT, nothing, nothing).time == PR.Infinity()
        @test PR.OptimalControlProblem(sys, XI, XT, nothing, nothing; safe_set = XS).safe_set ==
              XS
    end

    @testset "SafetyProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-3.0, -3.0), high = SVector(3.0, 3.0))
        sys = DummySystem(X)

        XI = LazySets.Hyperrectangle(; low = SVector(-0.5, -0.5), high = SVector(0.5, 0.5))
        XS = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        T = 7

        p = PR.SafetyProblem(sys, XI, XS, T)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.safe_set == XS
        @test p.time == T
    end

    @testset "CoSafeLTLProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
        sys = DummySystem(X)

        XI =
            LazySets.Hyperrectangle(; low = SVector(-0.9, -0.9), high = SVector(-0.8, -0.8))

        # SPEC can be anything; use a Symbol or String as stand-in
        spec = :F_target

        # Label payload can be sets (concrete) or abstract states;
        # here use sets to match your recipe expectation.
        lab = Dict{Symbol, Any}(
            :goal => LazySets.Hyperrectangle(;
                low = SVector(0.4, 0.4),
                high = SVector(0.6, 0.6),
            ),
            :avoid => LazySets.Hyperrectangle(;
                low = SVector(-0.2, -0.2),
                high = SVector(0.2, 0.2),
            ),
        )

        ap_sem = Dict{Symbol, Any}(:goal => MP.INNER, :avoid => MP.OUTER)

        p = PR.CoSafeLTLProblem(sys, XI, spec, lab, ap_sem)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.spec == spec
        @test p.labeling == lab
        @test p.ap_semantics == ap_sem
    end

    @testset "Problem hierarchy and interface" begin
        X = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        sys = DummySystem(X)
        XI =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(-0.5, -0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(1.0, 1.0))

        p = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 5)
        @test p isa PR.ControlProblem
        @test !(p isa PR.AbstractionProblem)

        # reach-avoid rounds the horizon down; the others round up
        @test PR.horizon_round_up(p) == false
        XS = LazySets.Hyperrectangle(; low = SVector(-1.5, -1.5), high = SVector(1.5, 1.5))
        @test PR.horizon_round_up(PR.SafetyProblem(sys, XI, XS, 5)) == true

        # remake swaps fields, re-infers type parameters, copies the rest
        XT2 = LazySets.Hyperrectangle(; low = SVector(0.6, 0.6), high = SVector(1.1, 1.1))
        p2 = PR.remake(p; target_set = XT2, time = 9)
        @test p2 isa PR.OptimalControlProblem
        @test p2.target_set == XT2
        @test p2.time == 9
        @test p2.initial_set === p.initial_set
        @test p2.system === p.system

        # `remake` walks the fields generically, so it reaches `safe_set` too.
        XS2 = LazySets.Hyperrectangle(; low = SVector(-1.8, -1.8), high = SVector(1.8, 1.8))
        @test PR.remake(p; safe_set = XS2).safe_set == XS2
        @test PR.remake(p; time = 3).safe_set === nothing
    end

    @testset "Horizon: Infinity and discretize_time" begin
        @test !isfinite(PR.Infinity())
        @test isinf(PR.Infinity())

        # round_up = true → ceil ("for at least T"); false → floor ("within at most T")
        @test PR.discretize_time(1.0, 0.3; round_up = true) == 4
        @test PR.discretize_time(1.0, 0.3; round_up = false) == 3
        @test PR.discretize_time(0.9, 0.3) == 3            # exact; default rounds up
        @test PR.discretize_time(PR.Infinity(), 0.3) === PR.Infinity()  # stays infinite
    end

    @testset "Default (infinite) horizon constructors" begin
        X = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        sys = DummySystem(X)
        XI =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(-0.5, -0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(1.0, 1.0))
        XS = LazySets.Hyperrectangle(; low = SVector(-1.5, -1.5), high = SVector(1.5, 1.5))

        @test PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0).time ===
              PR.Infinity()
        @test PR.SafetyProblem(sys, XI, XS).time === PR.Infinity()
        @test PR.ReachAndStayProblem(sys, XI, XT, XS).time === PR.Infinity()
    end

    @testset "ReachAndStayProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        sys = DummySystem(X)
        XI =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(-0.5, -0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(1.0, 1.0))
        XS = LazySets.Hyperrectangle(; low = SVector(-1.5, -1.5), high = SVector(1.5, 1.5))

        p = PR.ReachAndStayProblem(sys, XI, XT, XS, 6)
        @test p isa PR.ControlProblem
        @test p.system === sys
        @test p.initial_set == XI
        @test p.target_set == XT
        @test p.safe_set == XS
        @test p.time == 6
        @test PR.horizon_round_up(p) == true   # reach-and-stay is a "for at least T" spec
    end

    @testset "BisimulationQuotientProblem fields" begin
        X = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
        sys = DummySystem(X)
        state_set =
            LazySets.Hyperrectangle(; low = SVector(-0.5, -0.5), high = SVector(0.5, 0.5))
        regions = [
            LazySets.Hyperrectangle(;
                low = SVector(-0.4, -0.4),
                high = SVector(-0.1, -0.1),
            ),
            LazySets.Hyperrectangle(; low = SVector(0.1, 0.1), high = SVector(0.4, 0.4)),
        ]

        p = PR.BisimulationQuotientProblem(sys, state_set, regions)
        @test p isa PR.ProblemType
        @test p isa PR.AbstractionProblem
        @test p.system === sys
        @test p.state_set == state_set
        @test p.observation_regions === regions
    end

    @testset "trajectory_success semantics" begin
        XI =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(-0.5, -0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5, 0.5), high = SVector(1.0, 1.0))
        XS = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
        sys = DummySystem(XS)

        in_XI = SVector(-0.75, -0.75)
        mid = SVector(0.0, 0.0)
        in_XT = SVector(0.75, 0.75)
        unsafe = SVector(3.0, 3.0)

        traj(pts) = ST.Trajectory(pts)
        empty_traj = ST.Trajectory(SVector{2, Float64}[])

        @testset "reach-avoid (OptimalControlProblem)" begin
            p = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 5)
            @test PR.trajectory_success(p, traj([in_XI, mid, in_XT]))  # from XI, reaches XT
            @test !PR.trajectory_success(p, traj([in_XI, mid]))        # never reaches XT
            @test !PR.trajectory_success(p, traj([mid, in_XT]))        # does not start in XI
            @test !PR.trajectory_success(p, empty_traj)
        end

        @testset "reach-avoid with a safe set" begin
            # A corridor that excludes `mid`, so the middle of the direct route is off limits.
            corridor = LazySets.Hyperrectangle(;
                low = SVector(-1.0, -1.0),
                high = SVector(-0.25, -0.25),
            )
            keep_out = PR.OptimalControlProblem(
                sys,
                XI,
                XT,
                nothing,
                nothing,
                5;
                safe_set = corridor,
            )
            @test !PR.trajectory_success(keep_out, traj([in_XI, mid, in_XT]))

            wide = LazySets.Hyperrectangle(;
                low = SVector(-2.0, -2.0),
                high = SVector(2.0, 2.0),
            )
            ok = PR.OptimalControlProblem(sys, XI, XT, nothing, nothing, 5; safe_set = wide)
            @test PR.trajectory_success(ok, traj([in_XI, mid, in_XT]))
            @test !PR.trajectory_success(ok, traj([in_XI, unsafe, in_XT]))  # leaves the safe set

            # `safe U target`: what happens after the target is reached is not constrained.
            @test PR.trajectory_success(ok, traj([in_XI, mid, in_XT, unsafe]))
        end

        @testset "safety (SafetyProblem)" begin
            XIs = LazySets.Hyperrectangle(;
                low = SVector(-0.5, -0.5),
                high = SVector(0.5, 0.5),
            )
            p = PR.SafetyProblem(sys, XIs, XS, 5)
            @test PR.trajectory_success(p, traj([mid, in_XT, in_XI]))  # from XI, stays in XS
            @test !PR.trajectory_success(p, traj([mid, unsafe]))       # leaves XS
            @test !PR.trajectory_success(p, traj([in_XT, mid]))        # in_XT ∉ XIs
            @test !PR.trajectory_success(p, empty_traj)
        end

        @testset "reach-and-stay (ReachAndStayProblem)" begin
            p = PR.ReachAndStayProblem(sys, XI, XT, XS, 5)
            @test PR.trajectory_success(p, traj([in_XI, mid, in_XT]))    # ends in XT, stays safe
            @test !PR.trajectory_success(p, traj([in_XI, in_XT, mid]))   # ends outside XT
            @test !PR.trajectory_success(p, traj([in_XI, unsafe, in_XT]))  # leaves XS
            @test !PR.trajectory_success(p, empty_traj)
        end

        @testset "reach-and-stay with stay_on_first_entry" begin
            loose = PR.ReachAndStayProblem(sys, XI, XT, XS, 5)
            strict = PR.ReachAndStayProblem(sys, XI, XT, XS, 5; stay_on_first_entry = true)
            @test !loose.stay_on_first_entry
            @test strict.stay_on_first_entry
            # The keyword survives the 4-argument (infinite-horizon) constructor too.
            @test PR.ReachAndStayProblem(sys, XI, XT, XS; stay_on_first_entry = true).stay_on_first_entry

            # The case that separates the two readings: enter the target, leave, come back.
            # ◇□ accepts it — some suffix is in the target. "Stay from first entry" does not.
            reentry = traj([in_XI, in_XT, mid, in_XT])
            @test PR.trajectory_success(loose, reentry)
            @test !PR.trajectory_success(strict, reentry)

            # A run that arrives and holds satisfies both.
            settles = traj([in_XI, mid, in_XT, in_XT])
            @test PR.trajectory_success(loose, settles)
            @test PR.trajectory_success(strict, settles)

            # Never reaching the target fails the strict reading for want of a first entry.
            @test !PR.trajectory_success(strict, traj([in_XI, mid, mid]))
        end

        @testset "co-safe LTL (placeholder returns false)" begin
            lab = Dict{Symbol, Any}(:ap => XT)
            ap_sem = Dict{Symbol, Any}(:ap => MP.INNER)
            p = PR.CoSafeLTLProblem(sys, XI, :spec, lab, ap_sem)
            @test !PR.trajectory_success(p, traj([in_XI, in_XT]))       # not yet implemented
            @test !PR.trajectory_success(p, empty_traj)
        end
    end

    @testset "discretize_problem" begin
        sys = single_integrator()   # a real continuous-time plant (ẋ = u)
        XI = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(-0.5))
        XT = LazySets.Hyperrectangle(; low = SVector(0.5), high = SVector(1.0))

        # reach-avoid rounds the horizon down: floor(1.0 / 0.3) = 3
        p = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 1.0)
        pd = PR.discretize_problem(p, 0.3)
        @test pd isa PR.OptimalControlProblem
        @test pd.time == 3
        @test pd.initial_set === XI       # non-system fields copied verbatim
        @test pd.system !== sys           # the system was time-discretized
        @test pd.system isa MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem

        # safety rounds up: ceil(1.0 / 0.3) = 4
        ps = PR.SafetyProblem(
            sys,
            XI,
            LazySets.Hyperrectangle(; low = SVector(-1.5), high = SVector(1.5)),
            1.0,
        )
        @test PR.discretize_problem(ps, 0.3).time == 4

        # abstraction problems have no generic discretization → the stub errors
        alt = PR.AlternatingSimulationProblem(sys, XI)
        @test_throws ErrorException PR.discretize_problem(alt, 0.3)
    end

    @testset "Specifications" begin
        set1 = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(1.0, 1.0))
        set2 =
            LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(0.0, 0.0))

        # StateSpec: default INNER, explicit mode, and time-agnostic membership
        s1 = PR.StateSpec(set1)
        @test s1 isa PR.AbstractSpecification
        @test s1.incl_mode === UT.INNER
        @test PR.StateSpec(set1, UT.OUTER).incl_mode === UT.OUTER
        @test PR.satisfies(s1, SVector(0.5, 0.5))
        @test !PR.satisfies(s1, SVector(2.0, 2.0))
        @test PR.satisfies(s1, SVector(0.5, 0.5), 99.0)   # StateSpec ignores time

        # TimedSpec: base holds AND t ∈ [tmin, tmax]
        ts = PR.TimedSpec(s1, 1.0, 3.0)
        @test PR.satisfies(ts, SVector(0.5, 0.5), 2.0)
        @test !PR.satisfies(ts, SVector(0.5, 0.5), 5.0)   # outside the time window
        @test !PR.satisfies(ts, SVector(2.0, 2.0), 2.0)   # outside the set

        # HybridSpec, time-free mode dispatch (per-mode StateSpec)
        hs = PR.HybridSpec(Dict(1 => s1, 2 => PR.StateSpec(set2)))
        @test PR.satisfies(hs, SVector(0.5, 0.5), 1)
        @test PR.satisfies(hs, SVector(-0.5, -0.5), 2)
        @test !PR.satisfies(hs, SVector(0.5, 0.5), 2)     # not in mode 2's set
        @test !PR.satisfies(hs, SVector(0.5, 0.5), 99)    # absent mode → false

        # hybrid_reach_spec: parallel (state, time, mode) → mode-indexed timed spec
        H = PR.hybrid_reach_spec(
            [set1, set2],
            [
                LazySets.Hyperrectangle(; low = SVector(0.0), high = SVector(2.0)),
                LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(4.0)),
            ],
            [1, 2],
        )
        @test H isa PR.HybridSpec
        @test PR.satisfies(H, SVector(0.5, 0.5), 1.0, 1)    # mode 1, in set, t ∈ [0, 2]
        @test !PR.satisfies(H, SVector(0.5, 0.5), 3.0, 1)   # t = 3 ∉ [0, 2]
        @test PR.satisfies(H, SVector(-0.5, -0.5), 2.0, 2)  # mode 2, in set, t ∈ [1, 4]
        @test !PR.satisfies(H, SVector(0.5, 0.5), 1.0, 9)   # absent mode → false

        # validation: parallel lengths must match and modes must be unique
        @test_throws AssertionError PR.hybrid_reach_spec(
            [set1, set2],
            [LazySets.Hyperrectangle(; low = SVector(0.0), high = SVector(2.0))],
            [1, 2],
        )
        @test_throws AssertionError PR.hybrid_reach_spec(
            [set1, set2],
            [
                LazySets.Hyperrectangle(; low = SVector(0.0), high = SVector(2.0)),
                LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(4.0)),
            ],
            [1, 1],
        )
    end

    # A problem recipe draws the sets that define the specification, one labelled series each.
    # The labels *are* the contract — they tell the reader which region is which — so they are
    # asserted rather than merely checking that nothing threw. Run through the full Plots
    # pipeline; see `test/utils/plotting.jl` for why that is the level that exercises a recipe.
    @testset "Plots recipes" begin
        X = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
        sys = DummySystem(X)
        XI =
            LazySets.Hyperrectangle(; low = SVector(-0.9, -0.9), high = SVector(-0.8, -0.8))
        XT = LazySets.Hyperrectangle(; low = SVector(0.8, 0.8), high = SVector(0.9, 0.9))
        XS = LazySets.Hyperrectangle(; low = SVector(-0.7, -0.7), high = SVector(0.7, 0.7))
        region =
            LazySets.Hyperrectangle(; low = SVector(-0.5, -0.5), high = SVector(0.5, 0.5))

        labels(p; kw...) = [s[:label] for s in plot(p; kw...).series_list]

        @test labels(PR.AlternatingSimulationProblem(sys, region)) == ["Domain", "Region"]

        # A reach-avoid problem only grows a "Safe set" series when it actually has one.
        p_opt = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 3)
        @test labels(p_opt) == ["Domain", "Initial set", "Target set"]
        p_avoid =
            PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 3; safe_set = XS)
        @test labels(p_avoid) == ["Domain", "Safe set", "Initial set", "Target set"]

        @test labels(PR.SafetyProblem(sys, XI, XS, 4)) ==
              ["Domain", "Safe set", "Initial set"]

        @test labels(PR.ReachAndStayProblem(sys, XI, XT, XS, PR.Infinity())) ==
              ["Domain", "Safe set", "Target set", "Initial set"]

        # One series per observation region, numbered; the enclosing region is optional.
        p_bisim = PR.BisimulationQuotientProblem(sys, XS, [XI, XT])
        @test labels(p_bisim) == ["Region", "O 1", "O 2"]
        @test labels(p_bisim; plot_region = false) == ["O 1", "O 2"]

        # A co-safe LTL problem draws one series per atomic proposition. The labelling is a
        # `Dict`, so only the set of names is well defined, not their order.
        lab = Dict{Symbol, Any}(
            :goal => LazySets.Hyperrectangle(;
                low = SVector(0.2, 0.2),
                high = SVector(0.3, 0.3),
            ),
            :hazard => LazySets.Hyperrectangle(;
                low = SVector(-0.3, -0.3),
                high = SVector(-0.2, -0.2),
            ),
        )
        ap_sem = Dict{Symbol, Any}(:goal => MP.INNER, :hazard => MP.OUTER)
        p_ltl = PR.CoSafeLTLProblem(sys, XI, :spec, lab, ap_sem)
        @test Set(labels(p_ltl)) == Set(["Domain", "Initial set", "goal", "hazard"])

        plt = plot(p_ltl; ap_colors = Dict{Symbol, Any}(:goal => :magenta))
        goal = plt.series_list[findfirst(s -> s[:label] == "goal", plt.series_list)]
        @test goal[:fillcolor] == Plots.plot_color(:magenta)
    end
end

println("End problem tests")

end # module TestMain
