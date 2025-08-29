module FlowShopScheduling1D

using StaticArrays, MathematicalSystems, HybridSystems
using LinearAlgebra

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

struct FlowShopResetMap <: MathematicalSystems.AbstractMap
    domain::UT.HyperRectangle
    x_init::Vector{Float64}
    t_min::Float64
end

MathematicalSystems.apply(reset::FlowShopResetMap, state::AbstractVector) =
    vcat(reset.x_init, max(reset.t_min, state[end]))
MathematicalSystems.stateset(reset::FlowShopResetMap) = reset.domain

# """
#     make_cost_function(mode_weights::Vector{Float64}, t_nexttask_starts::Vector{Float64}; switch_penalty=100.0, base_switch_cost=1.0)

# Constructs a custom cost function for the 1D flowshop.

# - `mode_weights`: weights per mode (e.g., [3.0, 2.0, ...])
# - `t_nexttask_starts`: start time of each task (e.g., [1.0, 7.0, ...])
# - `switch_penalty`: penalty coefficient for switching (default 100.0)
# - `base_switch_cost`: base cost when switching (default 1.0)
# """
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

# """
#     generate_system_and_problem()

#     Generate a 1D flowshop scheduling hybrid control problem with 5 sequential tasks.

#     - State: (x, t, k) where x ∈ ℝ¹ (system state), t ∈ ℝ (time), k ∈ {1,2,3,4,5} (mode/task index)
#     - Each mode/task has its own continuous dynamics, state/input/time constraints, and guard (acceptance region).
#     - Guards: Each guard is a rectangle in (x, t) defining the acceptance region for switching to the next task.
#     - Reset maps: When a guard is reached, the state is reset (x, t) → (x_init, max(t_min, t)).
#     - The automaton encodes the allowed sequence of tasks (1→2→3→4→5).
#     - The final target is a region in (x, t) for the last mode.
#     - The cost function is mode-dependent, penalizes input effort, and strongly penalizes switching before the end of the time window (to encourage waiting as long as possible before switching).

#     Guards (acceptance regions):
#         - Task 1: x ∈ [6,10], t ∈ [0,3]
#         - Task 2: x ∈ [8,12], t ∈ [1,5]
#         - Task 3: x ∈ [10,11], t ∈ [7,9]
#         - Task 4: x ∈ [7,10], t ∈ [8,11]
#         - Task 5 (target): x ∈ [8,10], t ∈ [10,13]

#     The problem is designed to test temporal logic, switching, and optimal control in a simple 1D setting.
# """
function generate_system_and_problem()

    # Discretization parameters (dx, du, dt) for each task
    discretization_parameters = [
        (0.25, 0.25, 0.2),
        (0.25, 0.5, 0.1),
        (0.1, 0.1, 0.2),
        (0.1, 0.1, 0.1),
        (0.2, 0.25, 0.1),
    ]

    # Create optimizer factories for each mode using UniformGridAbstraction
    optimizer_factory_list = [
        () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
        () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
        () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
        () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
        () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    ]

    # Create state and input grids for each mode
    state_grid_1 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[1][1]))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[1][2]))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[2][1]))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[2][2]))
    state_grid_3 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[3][1]))
    input_grid_3 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[3][2]))
    state_grid_4 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[4][1]))
    input_grid_4 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[4][2]))
    state_grid_5 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[5][1]))
    input_grid_5 = DO.GridFree(SVector(0.0), SVector(discretization_parameters[5][2]))

    # Create optimizer parameters dictionary
    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => discretization_parameters[1][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => discretization_parameters[2][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
        Dict(
            "state_grid" => state_grid_3,
            "input_grid" => input_grid_3,
            "time_step" => discretization_parameters[3][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.6),
        ),
        Dict(
            "state_grid" => state_grid_4,
            "input_grid" => input_grid_4,
            "time_step" => discretization_parameters[4][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.7),
        ),
        Dict(
            "state_grid" => state_grid_5,
            "input_grid" => input_grid_5,
            "time_step" => discretization_parameters[5][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    # Dynamics for each task
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
    task1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task1_dynamics,
        1,
        1,
        X1,
        U1,
    )
    task2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task2_dynamics,
        1,
        1,
        X2,
        U2,
    )
    task3_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task3_dynamics,
        1,
        1,
        X3,
        U3,
    )
    task4_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task4_dynamics,
        1,
        1,
        X4,
        U4,
    )
    task5_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task5_dynamics,
        1,
        1,
        X5,
        U5,
    )

    # Time systems for each task
    timewindow_task1 = UT.HyperRectangle([0.0], [3.0]);
    task_1_time_system =
        MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task1)
    timewindow_task2 = UT.HyperRectangle([1.0], [5.0]);
    task_2_time_system =
        MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task2)
    timewindow_task3 = UT.HyperRectangle([7.0], [9.0]);
    task_3_time_system =
        MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task3)
    timewindow_task4 = UT.HyperRectangle([8.0], [11.0]);
    task_4_time_system =
        MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task4)
    timewindow_task5 = UT.HyperRectangle([10.0], [13.0]);
    task_5_time_system =
        MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task5)

    # Mode systems for the automaton
    modes_systems = [
        SY.VectorContinuousSystem([task1_system, task_1_time_system]),
        SY.VectorContinuousSystem([task2_system, task_2_time_system]),
        SY.VectorContinuousSystem([task3_system, task_3_time_system]),
        SY.VectorContinuousSystem([task4_system, task_4_time_system]),
        SY.VectorContinuousSystem([task5_system, task_5_time_system]),
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
    automaton = HybridSystems.GraphAutomaton(5)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    HybridSystems.add_transition!(automaton, 2, 3, 2)
    HybridSystems.add_transition!(automaton, 3, 4, 3)
    HybridSystems.add_transition!(automaton, 4, 5, 4)

    switchings = [
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
    ]

    HybridSystem_automaton =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Initial state and target
    initial_state = ([-0.5], 0.0, 1)
    Xs_target = [UT.HyperRectangle([8.0], [10.0])]
    Ts_target = [timewindow_task5]
    Ns_target = [5]

    mode_weights = [3.0, 11.0, 1.5, 1.2, 2.5]
    t_nexttask_starts = [1.0, 7.0, 8.0, 10.0]
    cost_function = make_cost_function(mode_weights, t_nexttask_starts)

    problem_specs = AB.TimedHybridAbstraction.TimedHybridOptimalControlProblem(
        initial_state,
        Xs_target,
        Ts_target,
        Ns_target,
        cost_function,
    )

    return HybridSystem_automaton,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    problem_specs
end
end
