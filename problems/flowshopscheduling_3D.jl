module FlowShopScheduling3D

# This module simulates a flowshop scheduling problem with 3 tasks and 3D dynamics (x, y, z).
# Each phase has its own linear dynamics and time window. The problem will be improved in the future with more interesting and realistic dynamics.

using LinearAlgebra 
using StaticArrays 

using Test
using StaticArrays, MathematicalSystems, HybridSystems

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

function create_jacobian(A_matrix)
    L = MMatrix{3,3,Float64}(A_matrix)
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

"""
    This function generates a 3-task flowshop scheduling problem with 3D linear dynamics.
    Each task has its own state and input constraints, time window, and reset map.
    The final target is a region in (x, y, z) at a given time window.
    This is a template for more advanced flowshop scheduling problems with richer dynamics.
"""
function generate_system_and_problem()

    # ------- Discretization parameters ------
    #(dx,  du,  dt)
    discretization_parameters = [
        (0.5, 0.5, 0.2), 
        (0.5, 0.5, 0.2), 
        (0.5, 0.5, 0.2), 
    ]

    # ------- Dynamic matrices inspired by a simplified 3D drone -------
    A1 = SMatrix{3,3}(0.95, 0.01, 0.003,
                      0.0,  0.98, 0.01,
                      -0.3,  0.0,  0.99)
    A2 = SMatrix{3,3}(0.90, 0.05, -0.1,
                      -0.01,  0.92, 0.03,
                      0.01, 0.0,  0.93)
    A3 = SMatrix{3,3}(0.85, 0.10, 0.0,
                      0.0,  0.88, -0.05,
                      0.02, 0.0,  0.90)

    B1 = SMatrix{3,3}(1.0, 0.1, 0.2,
                      0.0, 0.8, 0.0,
                      0.0, 0.0, 1.0)
    B2 = SMatrix{3,3}(0.9, 0.2, 0.1,
                      -0.001, 0.7, 0.1,
                      0.0, 0.1, 1.2)


    # ------- Adapted growth bounds -------
    growth_bounds = SVector(
        create_jacobian(A1),
        create_jacobian(A2),
        create_jacobian(A3),
    )

    # ------- Hybrid system definition -------
    task1_dynamics(x, u) = A1 * x + B1 * u
    task2_dynamics(x, u) = A2 * x + B2 * u
    task3_dynamics(x, u) = A3 * x .+ u

    X1 = UT.HyperRectangle([-1.0, -1.0, 0.0], [5.0, 5.0, 6.0]);
    U1 = UT.HyperRectangle([-1.0, -1.0, -0.5], [4.0, 4.0, 5.5]);
    X2 = UT.HyperRectangle([-2.5, -2.5, 0.0], [2.5, 2.5, 3.5]);
    U2 = UT.HyperRectangle([-1.2, -1.2, -0.7], [9.2, 6.2, 5.7]);
    X3 = UT.HyperRectangle([-1.0, -1.0, 0.0], [5.0, 3.0, 4.0]);
    U3 = UT.HyperRectangle([-1.5, -1.5, -1.0], [7.5, 5.5, 4.0]);

    task1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task1_dynamics, 3, 3, X1, U1)
    task2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task2_dynamics, 3, 3, X2, U2)
    task3_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        task3_dynamics, 3, 3, X3, U3)

    timewindow_task1 = UT.HyperRectangle([0.0], [2.0]);
    task_1_time_system = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task1)
    timewindow_task2 = UT.HyperRectangle([1.5], [4.0]);
    task_2_time_system = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task2)
    timewindow_task3 = UT.HyperRectangle([5.0], [7.0]);
    task_3_time_system = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], timewindow_task3)

    modes_systems = [
        SY.VectorContinuousSystem([task1_system, task_1_time_system]),
        SY.VectorContinuousSystem([task2_system, task_2_time_system]),
        SY.VectorContinuousSystem([task3_system, task_3_time_system]),
    ]
    task1_target = UT.HyperRectangle([0.5, 1.5, 1.0, 0.0], [5.0, 5.0, 6.0, 2.0])
    task2_target = UT.HyperRectangle([1.0, 0.0, 0.0, 1.5], [2.5, 2.5, 3.5, 4.0])
    task3_target = UT.HyperRectangle([0.5, 0.5, 0.0, 5.0], [5.0, 3.0, 4.0, 7.0])

    # Resets 
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0, 0.0, 2.0], 2.0)
    t2_t3_reset_map = FlowShopResetMap(task2_target, [1.0, 1.0, 3.0], 4.0)

    reset_maps = [
        t1_t2_reset_map,
        t2_t3_reset_map,
    ]

    # Transition automaton
    automaton = HybridSystems.GraphAutomaton(3)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    HybridSystems.add_transition!(automaton, 2, 3, 2)

    switchings = [
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
    ]

    HybridSystem_automaton = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Initial state (3D)
    initial_state = ([0.0, 0.0, 1.0], 0.0, 1)
    Xs_target = [UT.HyperRectangle([0.5, 0.5, 0.0], [5.0, 3.0, 4.0])]
    Ts_target = [timewindow_task3]
    Ns_target = [3]
    cost_function = (x, u) -> 1.0
    problem_specs = AB.TemporalHybridSymbolicModelAbstraction.ProblemSpecs(
        initial_state,
        Xs_target,
        Ts_target,
        Ns_target,
        cost_function,
    )

    return HybridSystem_automaton, growth_bounds, discretization_parameters, problem_specs
end
end