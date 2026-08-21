module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "AdaptiveCruiseControl",
        "adaptive_cruise_control.jl",
    ),
)

const ACC = AdaptiveCruiseControl
const params = ACC.Params()

function numerical_jacobian(f, x, u; h = 1e-6)
    e1, e2 = SVector(h, 0.0), SVector(0.0, h)
    return hcat((f(x + e1, u) - f(x - e1, u)) / (2h), (f(x + e2, u) - f(x - e2, u)) / (2h))
end

@testset "AdaptiveCruiseControl model" begin
    f = ACC.dynamic(params)
    J = ACC.jacobian(params)

    # At the lead's speed the gap is stationary, and holding that speed costs exactly the
    # resistance force — well inside the friction limit, or the benchmark would be trivial.
    @test f(SVector(50.0, params.v_lead), SVector(0.0))[1] ≈ 0.0
    a_hold = ACC.rolling_resistance(params, params.v_lead)
    @test f(SVector(50.0, params.v_lead), SVector(a_hold))[2] ≈ 0.0 atol = 1e-12
    @test 0.0 < a_hold < params.a_max

    # `SMatrix` fills column by column, so a transposed off-diagonal type-checks and runs.
    for x in (SVector(10.0, 5.0), SVector(80.0, 27.0))
        u = SVector(0.3)
        @test Matrix(J(x, u)) ≈ numerical_jacobian(f, x, u) atol = 1e-6
    end
end

@testset "AdaptiveCruiseControl sets" begin
    S = ACC.safe_set(params)

    @test SVector(50.0, 20.0) ∈ S       # 50 m ≥ 1.8 ⋅ 20 = 36 m
    @test SVector(30.0, 20.0) ∉ S       # 30 m < 36 m: tailgating
    @test SVector(50.0, 20.0) ∈ ACC.state_set(params)

    # Every target is inside the safe set, so reaching one cannot already have violated the
    # specification it was reached under.
    @test ACC.cruise_set(params) ⊆ S
    @test ACC.follow_set(params) ⊆ S

    # The set-point is faster than the lead; that is what makes the constraint bind.
    @test params.v_desired > params.v_lead
end

@testset "AdaptiveCruiseControl analytical bracket" begin
    a_outer = ACC.max_deceleration(params, params.v_max)
    @test a_outer > params.a_max

    # Below v⋆ the headway can be held with no margin at all; the margin is what bends the
    # boundary above it.
    v_star = params.v_lead + params.τ_h * params.a_max
    @test ACC.min_safe_gap(params, v_star) ≈ params.τ_h * v_star
    @test ACC.min_safe_gap(params, params.v_lead) ≈ params.τ_h * params.v_lead
    @test ACC.min_safe_gap(params, params.v_max) > params.τ_h * params.v_max

    for v in range(params.v_min, params.v_max; length = 50)
        inner = ACC.min_safe_gap(params, v)
        outer = ACC.min_safe_gap(params, v; a_brake = a_outer)
        @test outer <= inner + 1e-12        # stronger braking cannot need a larger gap
        @test inner >= params.τ_h * v
    end
end

@testset "AdaptiveCruiseControl safety synthesis" begin
    problem = ACC.safety_problem(; params = params)

    # Deliberately coarse: what is under test is the model against its closed form, not the
    # resolution the benchmark would be run at.
    hx = SVector(2.0, 0.25)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0, 0.0), hx),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.5)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        ACC.jacobian_bound(params),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    XMapping = SY.get_state_mapping(abstract_system)
    invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    @test MP.get_n_state(invariant_set, XMapping) > 0

    # Soundness against ground truth. No cell of the computed invariant set may sit below the
    # gap that even the strongest available braking needs — that region is unsafe as a matter
    # of physics, whatever the abstraction says. The slack covers one cell in each axis,
    # expressed in gap units.
    a_outer = ACC.max_deceleration(params, params.v_max)
    tol = hx[1] + params.τ_h * hx[2]
    worst = -Inf
    for q in MP.enum_states(invariant_set, XMapping)
        cell = MP.get_elem_by_state(XMapping, q)
        needed = ACC.min_safe_gap(params, LazySets.high(cell, 2); a_brake = a_outer)
        worst = max(worst, needed - LazySets.low(cell, 1))
    end
    @test worst <= tol

    # A degenerate box is the single point; OUTER returns the cell it falls in.
    function in_invariant(z, v)
        states = MP.get_states_from_set(
            XMapping,
            LazySets.Hyperrectangle(; low = SVector(z, v), high = SVector(z, v)),
            MP.OUTER,
        )
        return any(q -> MP.contains_state(invariant_set, XMapping, q), states)
    end

    @test in_invariant(60.0, 20.0)      # 24 m of margin over the headway at 20 m/s
    @test !in_invariant(10.0, 25.0)     # 10 m at 25 m/s, where 45 m are required
    # The envelope binds from above too: the gap opens at 8.9 m/s and the ego cannot reach the
    # lead's speed before running out of state space.
    @test !in_invariant(95.0, 5.0)

    x0 = SVector(60.0, 20.0)
    @test controller_admissible(controller, x0)
    traj = ST.get_closed_loop_trajectory(
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
        controller,
        x0,
        100,
    )
    @test PR.trajectory_success(problem, traj)
end

end # module TestMain
