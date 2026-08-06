module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using MathOptInterface
const MOI = MathOptInterface
import MathematicalSystems as MS
using HybridSystems

# Fast unit tests for the HybridSystemAbstraction helpers and dispatch logic —
# the parts the (slow) end-to-end suite either never touches (channelled_trajectory)
# or only hits on the happy path (the error branches of _find_switch_transition /
# control_solver_for). No abstraction is built here, so these stay fast.

const HSA = AB.HybridSystemAbstraction

# A tiny two-mode hybrid system with a single 1 → 2 transition.
function _tiny_hybrid_system()
    X = UT.box([-1.0], [1.0])
    U = UT.box([-1.0], [1.0])
    f(x, u) = x
    mode = MS.ConstrainedBlackBoxControlContinuousSystem(f, 1, 1, X, U)
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    return HybridSystems.HybridSystem(
        automaton,
        [mode, mode],
        [ST.GuardedResetMap(UT.box([-1.0], [1.0]))],
        [HybridSystems.AutonomousSwitching()],
    )
end

function _dummy_ocp()
    XI = UT.box(SVector(-1.0), SVector(-0.5))
    XT = UT.box(SVector(0.5), SVector(1.0))
    return PR.OptimalControlProblem(
        single_integrator(),
        XI,
        XT,
        nothing,
        (a, u) -> 1.0,
        PR.Infinity(),
    )
end

@testset "Optimizer construction & reset!" begin
    opt = HSA.Optimizer()
    @test opt isa HSA.Optimizer{Float64}
    @test opt.concrete_controller === nothing
    @test opt.print_level == 1
    @test opt.solve_time_sec == 0.0
    @test OP.default_abstraction_solver(opt) isa HSA.OptimizerAlternatingSimulationProblem

    # reset! clears outputs and forwards to both sub-solvers when they are present
    opt.control_solver = HSA.OptimizerOptimalControlProblem()
    opt.abstraction_solver = HSA.OptimizerAlternatingSimulationProblem()
    opt.solve_time_sec = 5.0
    returned = HSA.reset!(opt)
    @test returned === opt
    @test opt.concrete_controller === nothing
    @test opt.solve_time_sec == 0.0
end

@testset "control_solver_for dispatch" begin
    XI = UT.box(SVector(-1.0), SVector(-0.5))
    XT = UT.box(SVector(0.5), SVector(1.0))
    XS = UT.box(SVector(-1.0), SVector(1.0))
    sys = single_integrator()

    ocp = PR.OptimalControlProblem(sys, XI, XT, nothing, (a, u) -> 1.0, PR.Infinity())
    safety = PR.SafetyProblem(sys, XI, XS, PR.Infinity())
    unsupported = PR.ReachAndStayProblem(sys, XI, XT, XS, PR.Infinity())

    @test HSA.control_solver_for(ocp) isa HSA.OptimizerOptimalControlProblem
    @test HSA.control_solver_for(safety) isa HSA.OptimizerSafetyProblem
    @test_throws ErrorException HSA.control_solver_for(unsupported)
end

@testset "sub-optimizer MOI basics" begin
    ocp_opt = HSA.OptimizerOptimalControlProblem()
    @test MOI.is_empty(ocp_opt)                 # empty until the concrete problem is set
    ocp_opt.concrete_problem = _dummy_ocp()
    @test !MOI.is_empty(ocp_opt)
    @test MOI.get(ocp_opt, MOI.SolveTimeSec()) == 0.0
    @test HSA.reset!(ocp_opt) === ocp_opt

    safe_opt = HSA.OptimizerSafetyProblem()
    @test MOI.is_empty(safe_opt)
    @test MOI.get(safe_opt, MOI.SolveTimeSec()) == 0.0
    @test HSA.reset!(safe_opt) === safe_opt

    as_opt = HSA.OptimizerAlternatingSimulationProblem()
    @test MOI.is_empty(as_opt)
    @test MOI.get(as_opt, MOI.SolveTimeSec()) == 0.0
    @test HSA.reset!(as_opt) === as_opt

    # _validate_model errors on a missing required field, and passes once it is set
    @test_throws ErrorException HSA._validate_model(
        as_opt,
        [:alternating_simulation_problem],
    )
    as_opt.alternating_simulation_problem =
        PR.AlternatingSimulationProblem(single_integrator(), nothing)
    @test HSA._validate_model(as_opt, [:alternating_simulation_problem]) === nothing
end

@testset "_find_switch_transition" begin
    hs = _tiny_hybrid_system()

    tr, target_mode = HSA._find_switch_transition(hs, "SWITCH 1 -> 2")
    @test target_mode == 2
    @test HybridSystems.source(hs.automaton, tr) == 1
    @test HybridSystems.target(hs.automaton, tr) == 2

    # unparsable label → error
    @test_throws ErrorException HSA._find_switch_transition(hs, "not a switch label")
    # well-formed label but there is no 2 → 1 transition → error
    @test_throws ErrorException HSA._find_switch_transition(hs, "SWITCH 2 -> 1")
end

@testset "channelled_trajectory" begin
    u_traj = [[0.5], "SWITCH 1 -> 2"]

    # timed augmented states (x, t, mode) → the times channel is populated
    aug_timed = [([0.0], 0.0, 1), ([0.1], 0.1, 1), ([0.2], 0.0, 2)]
    traj = HSA.channelled_trajectory(aug_timed, u_traj)
    @test ST.states(traj) == [[0.0], [0.1], [0.2]]
    @test ST.modes(traj) == [1, 1, 2]
    @test ST.times(traj) == [0.0, 0.1, 0.0]
    @test ST.inputs(traj) == u_traj

    # time-free augmented states (x, mode) → no times channel
    aug_free = [([0.0], 1), ([0.1], 1), ([0.2], 2)]
    traj2 = HSA.channelled_trajectory(aug_free, u_traj)
    @test ST.states(traj2) == [[0.0], [0.1], [0.2]]
    @test ST.modes(traj2) == [1, 1, 2]
    @test ST.times(traj2) === nothing
end

@testset "reached / safe predicates" begin
    # reach spec: mode 2, x ∈ [-1, 1], t ∈ [1, 2]
    reach_spec = PR.hybrid_reach_spec(
        [UT.box(SVector(-1.0), SVector(1.0))],
        [UT.box(SVector(1.0), SVector(2.0))],
        [2],
    )
    ocp = PR.OptimalControlProblem(
        single_integrator(),
        ([0.0], 0.0, 1),
        reach_spec,
        nothing,
        (a, u) -> 1.0,
        PR.Infinity(),
    )
    @test HSA.reached(ocp, ([0.5], 1.5, 2))
    @test !HSA.reached(ocp, ([0.5], 1.5, 1))    # wrong mode
    @test !HSA.reached(ocp, ([2.0], 1.5, 2))    # x outside
    @test !HSA.reached(ocp, ([0.5], 3.0, 2))    # t outside the window

    safe_spec = PR.hybrid_reach_spec(
        [UT.box(SVector(-3.0), SVector(3.0))],
        [UT.box(SVector(0.0), SVector(5.0))],
        [1];
        incl_mode = UT.OUTER,
    )
    safety = PR.SafetyProblem(single_integrator(), ([0.0], 0.0, 1), safe_spec, 10.0)
    @test HSA.safe(safety, ([0.0], 1.0, 1))
    @test !HSA.safe(safety, ([4.0], 1.0, 1))     # x outside the safe set
end

end # module
