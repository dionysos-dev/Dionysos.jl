module TestStateInputVariationUGA

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
import MathOptInterface as MOI

# A state-dependent input constraint, written on concrete cells and inputs, goes
# through the UniformGridAbstraction front-end: here a "speed-limit zone" where only
# slow commands are allowed. The 1-D integrator lies on an exact lattice
# (`tstep * du == dx`), so the abstraction is deterministic.
@testset "state-dependent constraint through the grid abstraction" begin
    dx, du, tstep = 0.25, 0.5, 0.5

    F_sys(x, u) = SVector(u[1])
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0), high = SVector(2.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        1,
        1,
        _X_,
        _U_,
    )
    problem = PR.OptimalControlProblem(
        concrete_system,
        LazySets.Hyperrectangle(; low = SVector(-1.25), high = SVector(-0.75)),
        LazySets.Hyperrectangle(; low = SVector(0.75), high = SVector(1.25)),
        nothing,
        nothing,
        PR.Infinity(),
    )

    # Fast commands are forbidden in any cell reaching into x > 0.
    in_zone(cell) = LazySets.high(cell)[1] > 0.0
    limit = OPDS.StateInputVariation(
        (cell, u_prev, u) -> !in_zone(cell) || abs(u[1]) <= du;
        target_input = SVector(0.0),
        initial_input = SVector(0.0),
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0), SVector(dx)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(du)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bounded_input_variation"), limit)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    reached(x) = x ∈ problem.target_set
    traj = ST.get_closed_loop_trajectory(
        system,
        controller,
        SVector(-1.0),
        40;
        stopping = reached,
    )
    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))
    @test reached(xs[end])
    # Fast outside the zone, slow inside it.
    @test any(u -> abs(u[1]) > du, us)
    for (x, u) in zip(xs, us)
        x[1] + dx / 2 > 0.0 && @test abs(u[1]) <= du
    end
end

end # module
