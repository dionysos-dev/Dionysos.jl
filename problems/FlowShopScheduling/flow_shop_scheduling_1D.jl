module FlowShopScheduling1D

using StaticArrays

import HybridSystems as HS
import MathematicalSystems as MS

using LinearAlgebra

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PB = DI.Problem

struct FlowShopResetMap <: MS.AbstractMap
    domain::UT.HyperRectangle
    x_init::Vector{Float64}
    t_min::Float64
end

MS.apply(reset::FlowShopResetMap, state::AbstractVector) =
    vcat(reset.x_init, max(reset.t_min, state[end]))
MS.stateset(reset::FlowShopResetMap) = reset.domain

# Generate a 1D flowshop scheduling hybrid control problem with 5 sequential tasks.
#  - State: (x, t, k) where x ∈ ℝ¹ (system state), t ∈ ℝ (time), k ∈ {1,2,3,4,5} (mode/task index)
#  - Each mode/task has its own continuous dynamics, state/input/time constraints, and guard (acceptance region).
#  - Guards: Each guard is a rectangle in (x, t) defining the acceptance region for switching to the next task.
#  - Reset maps: When a guard is reached, the state is reset (x, t) → (x_init, max(t_min, t)).
#  - The automaton encodes the allowed sequence of tasks (1→2→3→4→5).
#  - The final target is a region in (x, t) for the last mode.
#  - The cost function is mode-dependent, penalizes input effort, and strongly penalizes switching before the end of the time window (to encourage waiting as long as possible before switching).
# Guards (acceptance regions):
#  - Task 1: x ∈ [6,10], t ∈ [0,3]
#  - Task 2: x ∈ [8,12], t ∈ [1,5]
#  - Task 3: x ∈ [10,11], t ∈ [7,9]
#  - Task 4: x ∈ [7,10], t ∈ [8,11]
#  - Task 5 (target): x ∈ [8,10], t ∈ [10,13]
# The problem is designed to test temporal logic, switching, and optimal control in a simple 1D setting.
function system()
    # Dynamics for each task/mode
    task1_dynamics(x, u) = [0.5 * x[1] + u[1]]
    task2_dynamics(x, u) = [0.8 * x[1] + u[1]]
    task3_dynamics(x, u) = [0.6 * x[1] + 0.6 * u[1]]
    task4_dynamics(x, u) = [0.7 * x[1] + 0.7 * u[1]]
    task5_dynamics(x, u) = task2_dynamics(x, u)

    # State and input spaces
    X1 = UT.HyperRectangle([-1.0], [10.0]);
    U1 = UT.HyperRectangle([-1.5], [5.5])
    X2 = UT.HyperRectangle([-1.0], [12.0]);
    U2 = UT.HyperRectangle([-1.5], [5.5])
    X3 = UT.HyperRectangle([1.0], [11.0]);
    U3 = UT.HyperRectangle([-1.5], [6.5])
    X4 = UT.HyperRectangle([0.0], [10.0]);
    U4 = UT.HyperRectangle([-1.5], [6.5])
    X5 = UT.HyperRectangle([-2.0], [10.0]);
    U5 = UT.HyperRectangle([-1.5], [5.5])

    # Continuous systems for each task
    task1_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task1_dynamics, 1, 1, X1, U1)
    task2_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task2_dynamics, 1, 1, X2, U2)
    task3_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task3_dynamics, 1, 1, X3, U3)
    task4_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task4_dynamics, 1, 1, X4, U4)
    task5_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task5_dynamics, 1, 1, X5, U5)

    # Time systems for each task
    timewindow_task1 = UT.HyperRectangle([0.0], [3.0]);
    task_1_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task1)
    timewindow_task2 = UT.HyperRectangle([1.0], [5.0]);
    task_2_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task2)
    timewindow_task3 = UT.HyperRectangle([7.0], [9.0]);
    task_3_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task3)
    timewindow_task4 = UT.HyperRectangle([8.0], [11.0]);
    task_4_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task4)
    timewindow_task5 = UT.HyperRectangle([10.0], [13.0]);
    task_5_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task5)

    # Mode systems for the automaton
    modes_systems = [
        ST.VectorContinuousSystem([task1_system, task_1_time_system]),
        ST.VectorContinuousSystem([task2_system, task_2_time_system]),
        ST.VectorContinuousSystem([task3_system, task_3_time_system]),
        ST.VectorContinuousSystem([task4_system, task_4_time_system]),
        ST.VectorContinuousSystem([task5_system, task_5_time_system]),
    ]

    # Guards (acceptance regions) for each task
    task1_target = UT.HyperRectangle([6.0, 0.0], [10.0, 3.0])
    task2_target = UT.HyperRectangle([8.0, 1.0], [12.0, 5.0])
    task3_target = UT.HyperRectangle([10.0, 7.0], [11.0, 9.0])
    task4_target = UT.HyperRectangle([7.0, 8.0], [10.0, 11.0])
    task5_target = UT.HyperRectangle([8.0], [10.0])

    # Reset maps for each transition
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0], 1.0)
    t2_t3_reset_map = FlowShopResetMap(task2_target, [2.0], 7.0)
    t3_t4_reset_map = FlowShopResetMap(task3_target, [1.0], 8.0)
    t4_t5_reset_map = FlowShopResetMap(task4_target, [-1.0], 10.0)

    reset_maps = [t1_t2_reset_map, t2_t3_reset_map, t3_t4_reset_map, t4_t5_reset_map]

    # Automaton with transitions between tasks
    automaton = HS.GraphAutomaton(5)
    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 3, 2)
    HS.add_transition!(automaton, 3, 4, 3)
    HS.add_transition!(automaton, 4, 5, 4)

    switchings = [
        HS.AutonomousSwitching(),
        HS.AutonomousSwitching(),
        HS.AutonomousSwitching(),
        HS.AutonomousSwitching(),
    ]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

function jacobian_bounds()
    return [
        u -> SMatrix{1, 1}(0.5),
        u -> SMatrix{1, 1}(0.8),
        u -> SMatrix{1, 1}(0.6),
        u -> SMatrix{1, 1}(0.7),
        u -> SMatrix{1, 1}(0.8),
    ]
end

# Constructs a custom cost function for the 1D flowshop.
# - `mode_weights`: weights per mode (e.g., [3.0, 2.0, ...])
# - `t_nexttask_starts`: start time of each task (e.g., [1.0, 7.0, ...])
# - `switch_penalty`: penalty coefficient for switching (default 100.0)
# - `base_switch_cost`: base cost when switching (default 1.0)
function make_cost_function(
    mode_weights::Vector{Float64},
    t_nexttask_starts::Vector{Float64};
    switch_penalty = 100.0,
    base_switch_cost = 1.0,
)
    return function (aug_state, u)
        (x, t, k) = aug_state
        w = mode_weights[k]
        if isa(u, String) && occursin("SWITCH", u)
            t_next = t_nexttask_starts[k]
            idle_penalty = (1+max(t_next - t, 0))^2
            return w * (base_switch_cost + switch_penalty * idle_penalty)
        end
        input_cost = (isa(u, Number) ? abs(u) : norm(u))
        return w * (1.0 + input_cost^2)
    end
end

function problem()
    # Hybrid System
    hybrid_system = system()

    # Initial state and target set
    initial_state = ([-0.5], 0.0, 1)

    Xs_target = [UT.HyperRectangle([8.0], [10.0])]
    Ts_target = [UT.HyperRectangle([10.0], [13.0])]
    Ns_target = [5]
    target_set = (Xs_target, Ts_target, Ns_target)

    # Cost Function
    mode_weights = [3.0, 11.0, 1.5, 1.2, 2.5]
    t_nexttask_starts = [1.0, 7.0, 8.0, 10.0]
    transition_cost_function = make_cost_function(mode_weights, t_nexttask_starts)

    return PB.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        transition_cost_function,
        PB.Infinity(),
    )
end

end # module
