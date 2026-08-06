module FlowShopScheduling3D

using StaticArrays

import HybridSystems as HS
import MathematicalSystems as MS

using LinearAlgebra

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

# On task completion the state restarts at `x_init`, and the clock is held at `t_min` when the
# task finished early. Guard and reset live in the augmented `[x; t]` space.
FlowShopResetMap(domain, x_init, t_min) =
    ST.GuardedResetMap(domain, state -> vcat(x_init, max(t_min, state[end])))

function get_matrices()
    A1 = SMatrix{3, 3}(0.95, 0.01, 0.003, 0.0, 0.98, 0.01, -0.3, 0.0, 0.99)
    A2 = SMatrix{3, 3}(0.90, 0.05, -0.1, -0.01, 0.92, 0.03, 0.01, 0.0, 0.93)
    A3 = SMatrix{3, 3}(0.85, 0.10, 0.0, 0.0, 0.88, -0.05, 0.02, 0.0, 0.90)
    return A1, A2, A3
end

# This function generates a 3-task flowshop scheduling problem with 3D linear dynamics.
# Each task has its own state and input constraints, time window, and reset map.
# The final target is a region in (x, y, z) at a given time window.
# This is a template for more advanced flowshop scheduling problems with richer dynamics.
function system()
    # Dynamics for each task/mode
    A1, A2, A3 = get_matrices()

    B1 = SMatrix{3, 3}(1.0, 0.1, 0.2, 0.0, 0.8, 0.0, 0.0, 0.0, 1.0)
    B2 = SMatrix{3, 3}(0.9, 0.2, 0.1, -0.001, 0.7, 0.1, 0.0, 0.1, 1.2)

    task1_dynamics(x, u) = A1 * x + B1 * u
    task2_dynamics(x, u) = A2 * x + B2 * u
    task3_dynamics(x, u) = A3 * x .+ u

    # State and input spaces
    X1 = UT.box([-1.0, -1.0, 0.0], [5.0, 5.0, 6.0])
    U1 = UT.box([-1.0, -1.0, -0.5], [4.0, 4.0, 5.5])
    X2 = UT.box([-2.5, -2.5, 0.0], [2.5, 2.5, 3.5])
    U2 = UT.box([-1.2, -1.2, -0.7], [9.2, 6.2, 5.7])
    X3 = UT.box([-1.0, -1.0, 0.0], [5.0, 3.0, 4.0])
    U3 = UT.box([-1.5, -1.5, -1.0], [7.5, 5.5, 4.0])

    # Continuous systems for each task
    task1_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task1_dynamics, 3, 3, X1, U1)
    task2_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task2_dynamics, 3, 3, X2, U2)
    task3_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task3_dynamics, 3, 3, X3, U3)

    # Time systems for each task
    timewindow_task1 = UT.box([0.0], [2.0])
    task_1_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task1)
    timewindow_task2 = UT.box([1.5], [4.0])
    task_2_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task2)
    timewindow_task3 = UT.box([5.0], [7.0])
    task_3_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task3)

    # Mode systems for the automaton
    modes_systems = [
        ST.VectorContinuousSystem([task1_system, task_1_time_system]),
        ST.VectorContinuousSystem([task2_system, task_2_time_system]),
        ST.VectorContinuousSystem([task3_system, task_3_time_system]),
    ]

    # Guards (acceptance regions) for each task
    task1_target = UT.box([0.5, 1.5, 1.0, 0.0], [5.0, 5.0, 6.0, 2.0])
    task2_target = UT.box([1.0, 0.0, 0.0, 1.5], [2.5, 2.5, 3.5, 4.0])
    task3_target = UT.box([0.5, 0.5, 0.0, 5.0], [5.0, 3.0, 4.0, 7.0])

    # Reset maps for each transition. `t_min` advances the clock to the start of the
    # next task's time window (task 2: 1.5, task 3: 5.0) so the switched state lands
    # inside the target task's clock domain.
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0, 0.0, 2.0], 2.0)
    t2_t3_reset_map = FlowShopResetMap(task2_target, [1.0, 1.0, 3.0], 5.0)

    reset_maps = [t1_t2_reset_map, t2_t3_reset_map]

    # Automaton with transitions between tasks
    automaton = HS.GraphAutomaton(3)
    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 3, 2)

    switchings = [HS.AutonomousSwitching(), HS.AutonomousSwitching()]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

function create_jacobian(A_matrix)
    L = MMatrix{3, 3, Float64}(A_matrix)
    n = size(L, 1)
    for i in 1:n
        for j in 1:n
            if i != j
                L[i, j] = abs(L[i, j])
            end
        end
    end
    return L
end

function jacobian_bounds()
    A1, A2, A3 = get_matrices()
    return [u -> create_jacobian(A1), u -> create_jacobian(A2), u -> create_jacobian(A3)]
end

function problem()
    # Hybrid System
    hybrid_system = system()

    # Initial state and target set
    initial_state = ([0.0, 0.0, 1.0], 0.0, 1)

    target_set = PR.hybrid_reach_spec(
        [UT.box([0.5, 0.5, 0.0], [5.0, 3.0, 4.0])],
        [UT.box([5.0], [7.0])],
        [3],
    )

    # Cost Function
    transition_cost_function = (x, u) -> 1.0

    return PR.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        transition_cost_function,
    )
end

end # module
