module FlowShopScheduling2D

using StaticArrays

import HybridSystems as HS
import MathematicalSystems as MS

using LinearAlgebra

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

# Reset on a guarded switch: reinitialize the 2-D state to `x_init` and push the clock
# forward to at least `t_min` (the state passed in is the augmented `[x1, x2, t]`).
struct FlowShopResetMap <: MS.AbstractMap
    domain::UT.Box
    x_init::Vector{Float64}
    t_min::Float64
end

MS.apply(reset::FlowShopResetMap, state::AbstractVector) =
    vcat(reset.x_init, max(reset.t_min, state[end]))
MS.stateset(reset::FlowShopResetMap) = reset.domain

function get_matrices()
    A1 = SMatrix{2, 2}(0.95, 0.01, 0.0, 0.98)
    A2 = SMatrix{2, 2}(0.90, 0.05, -0.01, 0.92)
    A3 = SMatrix{2, 2}(0.85, 0.10, 0.0, 0.88)
    return A1, A2, A3
end

# 3-task flowshop with 2-D linear dynamics per task plus a per-task clock. Each task
# has its own state/input constraints and time window; a guarded switch resets the
# state and advances the clock, and the final target is a region in (x1, x2) within a
# time window. This is the 2-D analogue of `flow_shop_scheduling_3D.jl`, small enough
# to abstract and solve on a laptop.
function system()
    A1, A2, A3 = get_matrices()

    B1 = SMatrix{2, 2}(1.0, 0.1, 0.0, 0.8)
    B2 = SMatrix{2, 2}(0.9, 0.2, -0.001, 0.7)

    task1_dynamics(x, u) = A1 * x + B1 * u
    task2_dynamics(x, u) = A2 * x + B2 * u
    task3_dynamics(x, u) = A3 * x .+ u

    # State and input spaces (2-D)
    X1 = UT.box([-1.0, -1.0], [5.0, 5.0]);
    U1 = UT.box([-1.0, -1.0], [4.0, 4.0]);
    X2 = UT.box([-2.5, -2.5], [2.5, 2.5]);
    U2 = UT.box([-1.2, -1.2], [9.2, 6.2]);
    X3 = UT.box([-1.0, -1.0], [5.0, 3.0]);
    U3 = UT.box([-1.5, -1.5], [7.5, 5.5]);

    task1_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task1_dynamics, 2, 2, X1, U1)
    task2_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task2_dynamics, 2, 2, X2, U2)
    task3_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task3_dynamics, 2, 2, X3, U3)

    # Time systems (clocks) for each task
    task_1_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], UT.box([0.0], [2.0]))
    task_2_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], UT.box([1.5], [4.0]))
    task_3_time_system = MS.ConstrainedLinearContinuousSystem([1.0;;], UT.box([5.0], [7.0]))

    modes_systems = [
        ST.VectorContinuousSystem([task1_system, task_1_time_system]),
        ST.VectorContinuousSystem([task2_system, task_2_time_system]),
        ST.VectorContinuousSystem([task3_system, task_3_time_system]),
    ]

    # Guards (acceptance regions) over [x1, x2, t]
    task1_target = UT.box([0.5, 1.5, 0.0], [5.0, 5.0, 2.0])
    task2_target = UT.box([1.0, 0.0, 1.5], [2.5, 2.5, 4.0])

    # Reset maps for each transition. The reset's `t_min` advances the clock to the
    # start of the next task's time window (task 2: 1.5, task 3: 5.0), so the switched
    # state lands inside the target task's clock domain.
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0, 0.0], 2.0)
    t2_t3_reset_map = FlowShopResetMap(task2_target, [1.0, 1.0], 5.0)
    reset_maps = [t1_t2_reset_map, t2_t3_reset_map]

    automaton = HS.GraphAutomaton(3)
    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 3, 2)
    switchings = [HS.AutonomousSwitching(), HS.AutonomousSwitching()]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

function create_jacobian(A_matrix)
    L = MMatrix{2, 2, Float64}(A_matrix)
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
    hybrid_system = system()

    initial_state = ([0.0, 0.0], 0.0, 1)   # (x1, x2, t, mode) -> ([x1,x2], t, mode)

    target_set =
        PR.hybrid_reach_spec([UT.box([0.5, 0.5], [5.0, 3.0])], [UT.box([5.0], [7.0])], [3])

    transition_cost_function = (aug_state, u) -> 1.0

    return PR.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        transition_cost_function,
    )
end

end # module
