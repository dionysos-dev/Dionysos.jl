module TestBipedRobot4DMultispeed

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

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
const OBSTACLE =
    LazySets.VPolygon([SVector(-0.39, 0.0), SVector(-0.35, 0.02), SVector(-0.31, 0.0)])
const X_BOX = LazySets.Hyperrectangle(;
    low = SVector(-0.6, -0.6, -0.6, -0.6),
    high = SVector(0.6, 0.6, 0.6, 0.6),
)

@testset "multi-speed steps with swept-cell validation" begin
    # speed_levels = 2: inputs up to two cells of displacement per axis per
    # step. Sound only with the swept filter, which rejects any transition
    # whose inter-sample segment crosses a removed cell.
    disc2 = RP.default_discretization(;
        dx = 0.1,
        tstep = 0.1,
        speed_levels = 2,
        swept_transitions = true,
    )
    removed = RP.infeasible_cells(X_BOX, disc2.state_grid, OBSTACLE, GEOM)
    domain = RP.RobotDomainConfig(;
        x_lb = SVector{4}(LazySets.low(X_BOX)),
        x_ub = SVector{4}(LazySets.high(X_BOX)),
        u_lb = SVector(-disc2.u_max, -disc2.u_max, -disc2.u_max, -disc2.u_max),
        u_ub = SVector(disc2.u_max, disc2.u_max, disc2.u_max, disc2.u_max),
    )
    sys = RP.system(;
        tstep = disc2.tstep,
        domain = domain,
        state_grid = disc2.state_grid,
        removed_cells = removed,
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(sys, nothing),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), disc2.state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), disc2.input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_input_filter"),
        MP.swept_input_filter(disc2.state_grid, disc2.tstep, removed),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("intersample_checked"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    @test SY.get_n_input(abstract_system) == 5^4

    # Independent soundness audit of the built transitions: no swept segment
    # crosses a removed cell, and multi-cell steps do occur.
    h = MP.get_h(disc2.state_grid)
    n_multi = 0
    n_violations = 0
    for (target, source, _) in SY.enum_transitions(abstract_system)
        xs_ = SY.get_concrete_state(abstract_system, source)
        xt_ = SY.get_concrete_state(abstract_system, target)
        # One-cell steps are sound by the base theorem (both endpoint cells are
        # kept); only multi-cell sweeps need the segment audit.
        maximum(abs.(xt_ - xs_) ./ h) > 1 + 1e-9 || continue
        n_multi += 1
        for pos in MP.cells_on_segment(disc2.state_grid, xs_, xt_)
            pos in removed && (n_violations += 1)
        end
    end
    @test n_multi > 0
    @test n_violations == 0

    x0 = SVector(0.2, 0.0, -0.2, 0.0)
    step_pb = RP.step_problem(sys, disc2.state_grid; x0 = x0)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

    # Slew-rate limited to one notch (du) per step: the controller must still
    # ramp through the intermediate speed level to exploit the fast inputs.
    # Re-setting `concrete_problem` after the constraint deliberately exercises
    # the attribute-preserving problem swap (the constraint used to be lost).
    slew = OPDS.BoundedInputVariation(
        (u1, u2) -> maximum(abs.(u1 - u2)),
        disc2.du;
        target_input = SVector(0.0, 0.0, 0.0, 0.0),
        initial_input = SVector(0.0, 0.0, 0.0, 0.0),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bounded_input_variation"), slew)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    dt_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    reached(x) = x ∈ step_pb.target_set
    traj = ST.get_closed_loop_trajectory(dt_system, controller, x0, 100; stopping = reached)
    txs = collect(ST.states(traj))
    tus = collect(ST.inputs(traj))
    @test reached(txs[end])
    # The richer alphabet is actually used, one notch at a time, rest to rest.
    @test maximum(maximum(abs.(u)) for u in tus) > disc2.du + 1e-9
    @test maximum(maximum(abs.(tus[k + 1] - tus[k])) for k in 1:(length(tus) - 1)) <=
          disc2.du + 1e-9
    @test maximum(abs.(tus[1])) <= disc2.du + 1e-9
    @test maximum(abs.(tus[end])) <= disc2.du + 1e-9
    # Inter-sample safety along the executed run, multi-cell steps included.
    n_bad = 0
    for k in 1:(length(txs) - 1)
        for pos in MP.cells_on_segment(disc2.state_grid, txs[k], txs[k + 1])
            pos in removed && (n_bad += 1)
        end
    end
    @test n_bad == 0
end

end # module
