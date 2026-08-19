module TestBipedRobot4DWalking

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI
import HybridSystems

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "BipedRobot",
        "4D_model",
        "robot_problem.jl",
    ),
)
import .RobotProblem as RP

const GEOM = RP.default_geometry()
const X_BOX = LazySets.Hyperrectangle(;
    low = SVector(-0.6, -0.6, -0.6, -0.6),
    high = SVector(0.6, 0.6, 0.6, 0.6),
)

@testset "leg swap is an involution and a grid automorphism" begin
    θ = SVector(0.11, -0.22, 0.33, -0.44)
    @test RP.leg_swap(θ) == SVector(0.33, -0.44, 0.11, -0.22)
    @test RP.leg_swap(RP.leg_swap(θ)) == θ

    # On a grid with one step for all axes the swap maps cell centres exactly onto
    # cell centres — this is what keeps the abstraction exact through the strike.
    grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), SVector(0.1, 0.1, 0.1, 0.1))
    for pos in ((1, -2, 3, 0), (-4, 4, -1, 2))
        c = SVector{4}(MP.get_coord_by_pos(grid, pos))
        swapped = RP.leg_swap(c)
        @test MP.get_coord_by_pos(grid, MP.get_pos_by_coord(grid, swapped)) ≈ swapped
        @test MP.get_pos_by_coord(grid, swapped) == (pos[3], pos[4], pos[1], pos[2])
    end
end

@testset "strike guard is a non-empty band of grounded, advanced cells" begin
    disc = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)
    band = RP.swing_foot_deviation_bound(GEOM, MP.get_h(disc.state_grid) ./ 2)
    guard = RP.strike_guard(
        X_BOX,
        disc.state_grid,
        GEOM;
        ground_band = band,
        min_advance = 0.05,
    )

    @test guard isa MP.CellUnion
    @test !isempty(guard)
    # Every guard cell really has its swing foot near the ground and ahead
    # (counted, so the record count does not grow with the guard).
    n_bad = 0
    for pos in guard
        foot = RP.swing_foot_position(
            SVector{4}(MP.get_coord_by_pos(disc.state_grid, pos)),
            GEOM,
            true,
        )
        (abs(foot[2]) <= band && foot[1] >= 0.05) || (n_bad += 1)
    end
    @test n_bad == 0
    # The guard is a proper subset: standing upright is not a strike.
    @test MP.get_pos_by_coord(disc.state_grid, SVector(0.0, 0.0, 0.0, 0.0)) ∉ guard
end

@testset "walking hybrid system: two modes, one shared abstraction" begin
    disc = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)
    domain = RP.RobotDomainConfig(;
        x_lb = SVector{4}(LazySets.low(X_BOX)),
        x_ub = SVector{4}(LazySets.high(X_BOX)),
        u_lb = SVector(-disc.u_max, -disc.u_max, -disc.u_max, -disc.u_max),
        u_ub = SVector(disc.u_max, disc.u_max, disc.u_max, disc.u_max),
    )
    sys = RP.system(;
        tstep = disc.tstep,
        domain = domain,
        state_grid = disc.state_grid,
        removed_cells = RP.infeasible_cells(X_BOX, disc.state_grid, nothing, GEOM),
    )
    guard = RP.strike_guard(X_BOX, disc.state_grid, GEOM)
    hs = RP.walking_hybrid_system(sys, guard)

    @test HybridSystems.nmodes(hs.automaton) == 2
    @test length(collect(HybridSystems.transitions(hs.automaton))) == 2
    # Both modes carry the very same system, which is what lets the two stance
    # phases share a single abstraction.
    @test HybridSystems.mode(hs, 1) === HybridSystems.mode(hs, 2)

    # The reset of each transition is the leg swap, applied to a guard cell.
    for t in HybridSystems.transitions(hs.automaton)
        reset = HybridSystems.resetmap(hs, t)
        θ = SVector(0.2, 0.0, -0.2, 0.1)
        @test MathematicalSystems.apply(reset, θ) == RP.leg_swap(θ)
        @test MathematicalSystems.stateset(reset) === guard
    end

    # The abstraction of the two modes is built once, not twice.
    kwargs = Dict(
        "state_grid" => disc.state_grid,
        "input_grid" => disc.input_grid,
        "time_step" => disc.tstep,
        "approx_mode" => AB.UniformGridAbstraction.CENTER_SIMULATION,
        "print_level" => 0,
    )
    calls = Ref(0)
    factory = function ()
        calls[] += 1
        return MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    end
    models = AB.HybridSystemAbstraction.build_mode_symbolic_models(
        hs,
        Function[factory, factory],
        [copy(kwargs), copy(kwargs)],
    )
    @test calls[] == 1
    @test models[2] === models[1]

    # Every guard cell resets to a cell that still exists in the domain, so no
    # strike edge is silently dropped.
    n_missing = 0
    for pos in guard
        swapped = RP.leg_swap(SVector{4}(MP.get_coord_by_pos(disc.state_grid, pos)))
        q = SY.get_abstract_state(models[1], swapped)
        (q === nothing || q <= 0) && (n_missing += 1)
    end
    @test n_missing == 0
end

@testset "gait recurrence: after a strike, the guard is reachable again" begin
    # Certifying walking is a *recurrence* property, not a safety one: staying in
    # the domain forever is trivially satisfied by standing still (u = 0 is an
    # input). What must hold is that the strike can be taken again and again.
    #
    # A finite computation certifies that infinite behaviour. Let `G` be the guard
    # and `Reach(G)` the states from which `G` is reachable. If every post-strike
    # state `π(G)` lies in `Reach(G)`, then by induction the robot can strike
    # forever — one reachability synthesis, no Büchi machinery.
    disc = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)
    w = RP.walking_problem(; geometry = GEOM, disc = disc, x_bar = 0.6, min_advance = 0.05)
    sol = RP.solve_walking(w)
    @test sol.success
    @test sol.controllable !== nothing

    # The recurrence certificate, computed by the model rather than restated here.
    n_post_strike, n_not_recurrent =
        RP.gait_recurrence(sol.abstract_system, w.guard, disc.state_grid, sol.controllable)
    @test n_post_strike > 0
    @test n_not_recurrent == 0

    # ------------------------------------------------------------------
    # Walk: chain the reachability controller with the strike and check that
    # several steps actually happen in closed loop.
    # ------------------------------------------------------------------
    X, guard, x0 = w.X, w.guard, w.x0
    aug_xs, us = RP.walk(w, sol, 120)

    strikes = count(u -> u isa AbstractString, us)
    @test strikes >= 3

    # The mode alternates at every strike and only there.
    modes = [aug[end] for aug in aug_xs]
    n_mode_changes = count(k -> modes[k + 1] != modes[k], 1:(length(modes) - 1))
    @test n_mode_changes == strikes

    # Every visited state stays in the carved domain, strikes included.
    removed = RP.infeasible_cells(X, disc.state_grid, nothing, GEOM)
    n_outside =
        count(aug -> MP.get_pos_by_coord(disc.state_grid, aug[1]) in removed, aug_xs)
    @test n_outside == 0

    # Each strike advances the robot: the swing foot lands ahead of the stance
    # foot, which is what makes the frame hop a step forward rather than in place.
    advances = Float64[]
    for (k, u) in enumerate(us)
        u isa AbstractString || continue
        foot = RP.swing_foot_position(SVector{4}(aug_xs[k][1]), GEOM, true)
        push!(advances, foot[1])
    end
    @test all(>=(0.05), advances)
end

end # module
