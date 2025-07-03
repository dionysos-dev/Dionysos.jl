# A time-augmented symbolic control framework around Dionysos

# TemporalHybridControl.jl

# This module extends the Dionysos symbolic control framework to support **time-augmented hybrid systems with task switching**.
# It allows the synthesis and execution of controllers across multiple temporally constrained tasks using symbolic abstractions.

module TemporalHybridControl

using StaticArrays, JuMP

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# -- Time mapping utilities --
struct TimeMapping
    tstep::Float64
end
time2int(tm::TimeMapping, t::Real) = floor(Int, t / tm.tstep)
int2time(tm::TimeMapping, idx::Int) = idx * tm.tstep
ceil_time2int(tm::TimeMapping, t::Real) = ceil(Int, t / tm.tstep)
floor_time2int(tm::TimeMapping, t::Real) = floor(Int, t / tm.tstep)

function add_task_transitions!(translist, task, k, abstract_system, tm)
    t_start = ceil_time2int(tm, task.a)
    t_end = floor_time2int(tm, task.b)

    for (target, source, abstract_input) in SY.get_transitions(abstract_system)
        for t in t_start:(t_end - 1)
            push!(translist, ((target, t + 1, k), (source, t, k), abstract_input))
        end
    end
end

function add_task_switches!(translist, tasks, abstract_system, tm)
    for k in 1:(length(tasks) - 1)
        task1 = tasks[k]
        task2 = tasks[k + 1]

        t1_abstract_target_set =
            SY.get_states_from_set(abstract_system, task1.target_set, DO.INNER), t_start =
                ceil_time2int(tm, max(task1.a, task2.b))
        t_end = floor_time2int(tm, task2.b)

        t2_initial_state = SY.get_abstract_state(abstract_system, task2.initial_state)
        t2_start = ceil_time2int(tm, task2.a)

        for t in t_start:t_end
            t_next = max(t, t2_start)
            for target in t1_abstract_target_set
                push!(translist, ((t2_initial_state, t_next, k + 1), (target, t, k), 0)) # 0 = input for task-switch
            end
        end
    end
end

function build_automaton(translist)
    function enum_augmented_states(transitions)
        states = Set{Any}()
        for (target, source, _) in transitions
            push!(states, target)
            push!(states, source)
        end
        return collect(states)
    end

    augmented_states = enum_augmented_states(translist)
    nstates = length(augmented_states)
    int2aug_state = [aug_state for aug_state in augmented_states]
    aug_state2int = Dict((aug_state, i) for (i, aug_state) in enumerate(augmented_states))

    autom = AutomatonList{S}(nx, nu)

    for (target, source, abstract_input) in translist
        target_int = aug_state2int[target]
        source_int = aug_state2int[source]
        SY.add_transition!(autom, source_int, target_int, abstract_input)
    end

    return int2aug_state, aug_state2int, autom
end

mutable struct TemporalHybridSymbolicModel{S1, A, N, T}
    symmodel::S1                              # Underlying symbolic model
    tm::T                                     # TimeMapping or similar type
    int2aug_state::Vector{NTuple{N, Int}}     # Integer → augmented state
    aug_state2int::Dict{NTuple{N, Int}, Int}  # Augmented state → integer
    autom::A                                  # Automaton representation (e.g., transition system)
end

function NewTemporalHybridSymbolicModel(symmodel, tasks, tstep::Float64)
    # -- Time discretization --
    tm = TimeMapping(tstep)

    translist = []
    # -- Build time-augmented transitions within task --
    for (k, task) in enumerate(tasks)
        add_task_transitions!(translist, task, k, symmodel, tm)
    end
    # -- Add task switches --
    add_task_switches!(translist, tasks, symmodel, tm)

    int2aug_state, aug_state2int, autom = build_automaton(translist)

    return TemporalHybridSymbolicModel(symmodel, tm, int2aug_state, aug_state2int, autom)
end

get_n_state(symmodel::TemporalHybridSymbolicModel) = length(symmodel.int2aug_state)
get_n_input(symmodel::TemporalHybridSymbolicModel) = SY.get_n_input(symmodel.symmodel)
enum_states(symmodel::TemporalHybridSymbolicModel) = 1:get_n_state(symmodel)
enum_inputs(symmodel::TemporalHybridSymbolicModel) = SY.enum_inputs(symmodel.symmodel)

function get_concrete_state(symmodel::TemporalHybridSymbolicModel, state)
    (q, t, k) = symmodel.int2aug_state[state]
    return (SY.get_concrete_state(symmodel.symmodel, q), int2time(symmodel.tm, t), k)
end

function get_concrete_input(symmodel::TemporalHybridSymbolicModel, input)
    return SY.get_concrete_input(symmodel.symmodel, input)
end
function get_abstract_state(symmodel::TemporalHybridSymbolicModel, aug_state)
    (x, t, k) = aug_state
    q = SY.get_abstract_state(symmodel.symmodel, x)
    t = floor_time2int(symmodel.tm, t)
    return (q, t, k)
end
function get_abstract_input(symmodel::TemporalHybridSymbolicModel, u)
    return SY.get_abstract_input(symmodel.symmodel, u)
end

function build_concrete_problem(tasks)
    init_state = tasks[1].initial_state
    concrete_initial_set = [(init_state, 0, 1)]
    concrete_target_set =
        (tasks[end].target_set, [tasks[end].a, tasks[end].b], length(tasks))

    return PR.OptimalControlProblem(
        tasks,
        concrete_initial_set,
        concrete_target_set,
        nothing,
        nothing,
        PB.Infinity(),
    )
end

function get_states_from_set(
    symmodel::TemporalHybridSymbolicModel,
    set,
    incl_mode::DO.INCL_MODE,
)
    (X, T, N) = set
    q_list = SY.get_states_from_set(symmodel.symmodel, X, incl_mode)
    abstract_time_interval =
        [ceil_time2int(symmodel.tm, T[1]), floor_time2int(symmodel.tm, T[2])]
    abstract_set = []
    for q in q_list
        for t in abstract_time_interval
            for k in N[1]:N[2]
                push!(abstract_set, get_abstract_state(symmodel, (q, t, k)))
            end
        end
    end
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    symmodel::TemporalHybridSymbolicModel,
)
    aug_state = concrete_problem.initial_set[1]
    abstract_initial_set = [get_abstract_state(symmodel, aug_state)]

    target_set = concrete_problem.target_set
    abstract_target_set = get_states_from_set(symmodel, target_set, DO.INNER)

    return PR.OptimalControlProblem(
        symmodel,
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost,       # TODO: Transform continuous cost into discrete abstraction
        concrete_problem.transitison_cost,  # TODO: Transform continuous cost into discrete abstraction
        concrete_problem.time,              # TODO: Translate continuous time into discrete steps
    )
end

function solve_abstract_problem(abstract_problem)
    abstract_controller,
    controllable_set_symbols,
    uncontrollable_set_symbols,
    value_fun_tab = OP.compute_largest_controllable_set(
        abstract_problem.abstract_system;
        target_set = abstract_problem.target_set,
        initial_set = abstract_problem.initial_set,
    )
    return abstract_controller
end

function solve_concrete_problem(symmodel::TemporalHybridSymbolicModel, abstract_controller)
    function concrete_controller(aug_state)
        (x, t, k) = aug_state
        abstract_aug_state = get_abstract_state(symmodel, aug_state)
        abstract_input =
            Dionysos.Utils.fix_and_eliminate_first(abstract_controller, abstract_aug_state)
        if abstract_input === 0
            return "Go to the next task"
        else
            return get_concrete_input(symmodel, abstract_input)
        end
    end
end

function solve(tasks, symmodel, tstep)
    temporalHybridSymmodel = NewTemporalHybridSymbolicModel(symmodel, tasks, tstep)
    # -- Build concrete problem --
    concrete_problem = build_concrete_problem(tasks)

    # -- Build abstract problem --
    abstract_problem = build_abstract_problem(concrete_problem, temporalHybridSymmodel)

    # -- Solve the abstract problem --
    abstract_controller = solve_abstract_problem(abstract_problem)

    # -- Solve the concrete problem using the abstract controller --
    concrete_controller =
        solve_concrete_problem(temporalHybridSymmodel, abstract_controller)

    return concrete_controller
end

function get_next_aug_state(tasks, continuous_time_system, aug_state, u, tstep)
    tm = TimeMapping(tstep)
    (x, t, k) = aug_state
    if u == "Go to the next task"
        return (
            tasks[k + 1].initial_state,
            max(t, ceil_time2int(tm, tasks[k + 1].a)),
            k + 1,
        )
    else
        next_x = ST.mapping(continuous_time_system, x, u, tstep)
        next_t = t + tstep
        return (next_x, next_t, k)
    end
end

function get_closed_loop_trajectory(
    tasks,
    continuous_time_system::MS.ConstrainedBlackBoxControlContinuousSystem,
    controller,
    aug_state_0,
    nstep;
    stopping = (x) -> false,
)
    aug_state_traj, u_traj = [aug_state_0], []
    aug_state = aug_state_0
    for _ in 1:nstep
        stopping(tasks, aug_state) && break
        u = controller(aug_state)
        u === nothing && break
        aug_state = get_next_aug_state(tasks, continuous_time_system, aug_state, u, tstep)
        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end
    return aug_state_traj, u_traj
end

function reached(tasks, aug_state)
    (x, t, k) = aug_state
    if k == length(tasks) && tasks[k].a <= t <= tasks[k].b && x ∈ tasks[k].target_set
        return true
    else
        return false
    end
end
end

function build_tasks()
    # ................................................
    # return tasks
end

function build_initial_abstraction(tasks, continuous_time_system, tstep)
    # ................................................
    # return symmodel
end

##### SCRIPTS 

# task.a => Float64
# task.b => Float64
# task.initial_state => SVector
# task.target_set => hyperrectangle
# tasks = [task1, task2, task3]  # Example tasks
# continuous_time_system => MS.ConstrainedBlackBoxControlContinuousSystem
# tstep => Float64

tasks = build_tasks()
symmodel = build_initial_abstraction(tasks, continuous_time_system, tstep)

concrete_controller = TemporalHybridControl.solve(tasks, symmodel, tstep)

aug_state_traj, u_traj = get_closed_loop_trajectory(
    tasks,
    concrete_controller,
    (tasks[1].initial_state, 0.0, 1);
    stopping = reached,
)
