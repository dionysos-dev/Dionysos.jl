module FlowShopScheduling1D

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
    domain::UT.HyperRectangle  # guard 
    x_init::Vector{Float64}    # reset point for X
    t_min::Float64             # t_min for the reset
end

MathematicalSystems.apply(reset::FlowShopResetMap, state::AbstractVector) =
    vcat(reset.x_init, max(reset.t_min, state[end]))
MathematicalSystems.stateset(reset::FlowShopResetMap) = reset.domain

""" 
    ICI il faudra expliquer le problème généré
"""
function generate_system_and_problem()

    # ------- Discretization parameters ------
    #(dx,  du,  dt)
    discretization_parameters = [
        (0.25, 0.25, 0.2), # task 1
        (0.25, 0.5, 0.1), # task 2
        (0.1, 0.1, 0.2),  # task 3
        (0.1, 0.1, 0.1),  # task 4
        (0.2, 0.25, 0.1),
    ] # task 5

    # ------- Define trivial growth bound ----
    growth_bounds = SVector(
        SMatrix{1, 1}(0.5),
        SMatrix{1, 1}(0.8),
        SMatrix{1, 1}(0.6),
        SMatrix{1, 1}(0.7),
        SMatrix{1, 1}(0.8),
    )

    # ------- Define the hybrid system -------
    # definition of the dynamics for each task
    task1_dynamics(x, u) = [0.5 * x[1] + u[1]]
    task2_dynamics(x, u) = [0.8 * x[1] + u[1]]
    task3_dynamics(x, u) = [0.6 * x[1] + 0.6 * u[1]]
    task4_dynamics(x, u) = [0.7 * x[1] + 0.7 * u[1]]
    task5_dynamics(x, u) = task2_dynamics(x, u)

    # define state_space and input_space for each task
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

    # creation of ConstrainedBlackBoxControlContinuousSystem for each task
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

    # create time system for each task
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
        SY.VectorContinuousSystem([task1_system, task_1_time_system]), # task 1
        SY.VectorContinuousSystem([task2_system, task_2_time_system]), # task 2
        SY.VectorContinuousSystem([task3_system, task_3_time_system]), # task 3
        SY.VectorContinuousSystem([task4_system, task_4_time_system]), # task 4
        SY.VectorContinuousSystem([task5_system, task_5_time_system]), # task 5
    ]

    # Define target for each task <-> guard of the mode (or global objective for the last task).
    task1_target = UT.HyperRectangle([6.0, 0.0], [10.0, 3.0]) # guard : [xmin_accept, tmin_window] x [xmax_accept, tmax_window]
    task2_target = UT.HyperRectangle([8.0, 1.0], [12.0, 5.0])
    task3_target = UT.HyperRectangle([5.0, 7.0], [11.0, 9.0])
    task4_target = UT.HyperRectangle([4.0, 8.0], [10.0, 11.0])
    task5_target = UT.HyperRectangle([4.0, 10.0], [10.0, 13.0])

    # Define the initial condition after finishing one task <-> the reset map associated with each mode.
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0], 3.0) # reset X to 0.0 and t to max(3.0, t) after finishing task 1 and before starting task 2
    t2_t3_reset_map = FlowShopResetMap(task2_target, [2.0], 5.0)
    t3_t4_reset_map = FlowShopResetMap(task3_target, [1.0], 9.0)
    t4_t5_reset_map = FlowShopResetMap(task4_target, [-1.0], 11.0)

    reset_maps = [
        t1_t2_reset_map, # reset map for transition from task 1 to task 2
        t2_t3_reset_map, # reset map for transition from task 2 to task 3
        t3_t4_reset_map, # reset map for transition from task 3 to task 4
        t4_t5_reset_map, # reset map for transition from task 4 to task 5
    ]

    # Define automaton with transitions between tasks
    automaton = HybridSystems.GraphAutomaton(5)
    HybridSystems.add_transition!(automaton, 1, 2, 1) # transition from task 1 to task 2
    HybridSystems.add_transition!(automaton, 2, 3, 2) # transition from task 2 to task 3
    HybridSystems.add_transition!(automaton, 3, 4, 3) # transition from task 3 to task 4
    HybridSystems.add_transition!(automaton, 4, 5, 4) # transition from task 4 to task 5

    switchings = [
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
    ]

    HybridSystem_automaton =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Define the initial state and input for the first task
    initial_state = ([-0.5], 0.0, 1) # initial state : task 1, x = -0.5, t = 0.0
    Xs_target = [task5_target] # target for the last task
    Ts_target = [timewindow_task5] # time window for the last task
    Ns_target = [5] # id of the mode target
    cost_function = (x, u) -> 1.0 # trivial cost function
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
