# Golden-output regression test.
#
# Pins *deterministic* results of the flagship UniformGridAbstraction pipeline so behaviour-preserving
# refactors are caught if they silently change the synthesized abstraction/controller.
# We pin only order-independent, randomness-free quantities (state/transition/controllable counts,
# success, reaches-target) — not the exact control input, which may be randomized.
#
# Regenerate the golden values after an *intentional* change:
#   DIONYSOS_SHOW_SIG=1 julia --project=test test/regression/uniform_grid_abstraction.jl
# then copy the printed numbers into the `@test` lines below.

module TestRegressionUniformGrid

using Test
using StaticArrays, MathematicalSystems
using Dionysos
using JuMP
import MathOptInterface as MOI

const DI = Dionysos
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/PathPlanning/path_planning.jl")

const SHOW_SIG = get(ENV, "DIONYSOS_SHOW_SIG", "0") == "1"

# Size a controllable set (an `AbstractStateSet`, e.g. `ExplicitIdSet`) via the mapping-aware API.
_set_size(s, mapping) = s === nothing ? 0 : MP.get_n_state(s, mapping)

@testset "regression: UniformGrid reachability (path planning, simple)" begin
    concrete_problem = PathPlanning.problem(; simple = true)

    state_grid = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2))
    input_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3))

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        PathPlanning.jacobian_bound(),
    )
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
    success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

    Xmap = SY.get_state_mapping(abstract_system)
    n_state = MP.get_n_state(Xmap)
    n_transitions = SY.get_n_transitions(abstract_system)
    n_controllable = _set_size(controllable_set, Xmap)

    # deterministic closed-loop reachability invariant
    reached(x) = (x ∈ concrete_problem.target_set)
    x0 = SVector(0.4, 0.4, 0.0)
    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        concrete_controller,
        x0,
        100;
        stopping = reached,
    )
    reaches_target = any(reached, collect(ST.enum_elems(traj.x)))

    if SHOW_SIG
        @info "UniformGrid reachability signature" n_state n_transitions n_controllable success reaches_target typeof(
            controllable_set,
        )
    end

    # ---- golden values (regenerate with DIONYSOS_SHOW_SIG=1) ----
    @test success == true
    @test reaches_target == true
    @test n_state == 26285
    @test n_transitions == 9705631
    @test n_controllable == 19703
end

end # module TestRegressionUniformGrid
