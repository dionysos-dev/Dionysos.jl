export TemporalHybridSymbolicModelAbstraction

module TemporalHybridSymbolicModelAbstraction

import HybridSystems
import HybridSystems: HybridSystem
import MathematicalSystems
import StaticArrays: SVector
import Dionysos

struct ProblemSpecs{F} # ajouter le type de problème
    initial_state::Tuple{AbstractVector{Float64}, Float64, Int} # aug_state initial ([x], t, mode_id)
    Xs_target::Vector{<:Dionysos.Utils.HyperRectangle} # must be changed to admit more complex target sets
    Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}
    Ns_target::Vector{Int}
    concret_cost_fun::F
end

# =================================
# CONTROLLER SYNTHESIS
# =================================

"""Builds the concrete control problem"""
function build_concrete_problem(Concret_specifications::ProblemSpecs)
    concrete_initial_set = Concret_specifications.initial_state
    concrete_target_set = (
        Concret_specifications.Xs_target,
        Concret_specifications.Ts_target,
        Concret_specifications.Ns_target,
    )
    concret_cost_fun = Concret_specifications.concret_cost_fun

    return Dionysos.Problem.OptimalControlProblem( # safety potentially
        Concret_specifications,
        concrete_initial_set,
        concrete_target_set,
        nothing,
        concret_cost_fun,
        Dionysos.Problem.Infinity(),
    )
end

"""Builds the abstract control problem"""
function build_abstract_problem(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel,
)
    aug_state = concrete_problem.initial_set
    abstract_initial_set =
        [Dionysos.Symbolic.TimedHybridAutomata.get_abstract_state(symmodel, aug_state)]
    target_set = concrete_problem.target_set
    abstract_target_set = Dionysos.Symbolic.TimedHybridAutomata.get_states_from_set(
        symmodel,
        target_set[1],
        target_set[2],
        target_set[3],
    )

    return Dionysos.Problem.OptimalControlProblem( # safety potentially
        symmodel,
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost,
        build_abstract_transition_cost(symmodel, concrete_problem.transition_cost),
        concrete_problem.time,
    )
end
"""Solves the abstract problem and computes the controllable set"""
function solve_abstract_problem(abstract_problem) # à voir si il faut changer 
    abstract_controller, controllable_set_symbols, _, value_per_node =
        Dionysos.Symbolic.compute_worst_case_cost_controller(
            abstract_problem.system.autom,
            abstract_problem.target_set;
            initial_set = abstract_problem.initial_set,
            sparse_input = false,
            cost_function = abstract_problem.transition_cost,
        )

    println("value init_state : ", value_per_node[abstract_problem.initial_set[1]])

    if ⊆(abstract_problem.initial_set, controllable_set_symbols)
        println("✅ Problem is solvable: initial set is controllable")
    else
        println("⚠️ Warning: initial set is only partially controllable")
    end

    return abstract_controller, controllable_set_symbols
end

"""Generates the concrete controller from the abstract controller"""
function solve_concrete_problem(
    symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel,
    abstract_controller,
)
    function concrete_controller(aug_state)
        (x, t, k) = aug_state
        abstract_aug_state =
            Dionysos.Symbolic.TimedHybridAutomata.get_abstract_state(symmodel, aug_state)

        if Dionysos.System.is_defined(abstract_controller, abstract_aug_state) == false
            println("⚠️ No action available for state $aug_state")
            return nothing
        end

        abstract_input =
            Dionysos.System.get_control(abstract_controller, abstract_aug_state)

        if Dionysos.Symbolic.TimedHybridAutomata.is_switching_input(
            symmodel.global_input_map,
            abstract_input,
        )
            transition_id = symmodel.global_input_map.global_to_switching[abstract_input]
            label = symmodel.global_input_map.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            return label
        else
            return Dionysos.Symbolic.TimedHybridAutomata.get_concrete_input(
                symmodel,
                abstract_input,
                k,
            )
        end
    end

    isdef = (aug_state) -> true # for the moment we assume the controller is defined for all augmented states
    return Dionysos.System.BlackBoxContinuousController(concrete_controller, isdef)
end

# =================================
# ABSTRACT COST FUNCTIONS
# =================================

function build_abstract_transition_cost(
    symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel,
    concrete_cost_fun,
)
    function abstract_transition_cost(state_int, input_int)
        (x, t, k) =
            Dionysos.Symbolic.TimedHybridAutomata.get_concrete_state(symmodel, state_int)
        aug_state = (x, t, k)
        if Dionysos.Symbolic.TimedHybridAutomata.is_switching_input(
            symmodel.global_input_map,
            input_int,
        )
            transition_id = symmodel.global_input_map.global_to_switching[input_int]
            label = symmodel.global_input_map.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            u = label
        else
            u = Dionysos.Symbolic.TimedHybridAutomata.get_concrete_input(
                symmodel,
                input_int,
                k,
            )
        end
        return concrete_cost_fun(aug_state, u)
    end
    return abstract_transition_cost
end

# =================================
# MAIN SOLVING FUNCTION
# =================================

"""
Solves the complete temporal hybrid control problem.

# Arguments
- `tasks`: List of tasks 
- `symmodel`: Base symbolic model
- `tstep`: Time discretization step

# Returns
Concrete controller for augmented states
"""
function solve(
    hs::HybridSystem,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    concret_specs::ProblemSpecs,
)
    hybrid_symmodel = Dionysos.Symbolic.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    concrete_problem = build_concrete_problem(concret_specs)

    abstract_problem = build_abstract_problem(concrete_problem, hybrid_symmodel)

    abstract_controller, _ = solve_abstract_problem(abstract_problem)

    concrete_controller = solve_concrete_problem(hybrid_symmodel, abstract_controller)

    return concrete_controller
end

# =================================
# CLOSED-LOOP SIMULATION
# =================================

"""Computes the next augmented state"""
function get_next_aug_state(hs::HybridSystem, aug_state, u, time_is_active, tstep, map_sys)
    (x, t, k) = aug_state

    if isa(u, AbstractString) && occursin("SWITCH", u)
        m = match(r"SWITCH (\d+) -> (\d+)", String(u))
        if m === nothing
            error("Unrecognized SWITCH label format: $u")
        end
        source_mode = parse(Int, m.captures[1])
        target_mode = parse(Int, m.captures[2])

        transitions = collect(HybridSystems.transitions(hs.automaton))
        transition_id = findfirst(
            tr ->
                HybridSystems.source(hs.automaton, tr) == source_mode &&
                HybridSystems.target(hs.automaton, tr) == target_mode,
            transitions,
        )
        if transition_id === nothing
            error("Transition not found for $u")
        end
        transition = transitions[transition_id]

        # Apply the reset map of the transition
        reset_map = HybridSystems.resetmap(hs, transition)
        augmented_source_state = vcat(x, t)
        reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
        next_x = reset_result[1:(end - 1)]
        next_t = reset_result[end]
        next_k = target_mode
        # Explicit rounding of time to 10 decimals to avoid error propagation
        next_t = round(next_t; digits = 10)
        return (next_x, next_t, next_k)
    else
        next_t = time_is_active ? t + tstep : 0.0
        # Explicit rounding of time to 10 decimals to avoid error propagation
        next_t = round(next_t; digits = 10)
        next_x = map_sys(x, u, tstep)
        return (next_x, next_t, k)
    end
end

"""
Generates a closed-loop trajectory.

# Arguments
- `tasks`: List of tasks
- `continuous_time_system`: Continuous-time dynamic system
- `controller`: Synthesized controller
- `aug_state_0`: Initial augmented state
- `nstep`: Maximum number of steps
- `stopping`: Stopping function

# Returns
Tuple (state_trajectory, control_trajectory)
"""
function get_closed_loop_trajectory(
    discretization_parameters::Vector{Tuple{Float64, Float64, Float64}},
    hs::HybridSystem,
    problem_specs::ProblemSpecs,
    controller,
    aug_state_0,
    nstep;
    stopping = (x) -> false,
)
    aug_state_traj, u_traj = [aug_state_0], []
    aug_state = aug_state_0

    nmodes = HybridSystems.nmodes(hs.automaton)
    dynamics = [HybridSystems.mode(hs, k).systems[1].f for k in 1:nmodes]
    maps_sys = [Dionysos.System.simulate_control_map(dynamics[i]) for i in 1:nmodes]
    tsteps = [discretization_parameters[i][3] for i in 1:nmodes]
    times_is_active =
        [([1.0;;]==HybridSystems.mode(hs, k).systems[2].A) ? true : false for k in 1:nmodes]

    for _ in 1:nstep
        stopping(problem_specs, aug_state) && break

        u = controller.f(aug_state)

        u === nothing && break
        (_, _, k) = aug_state
        aug_state =
            get_next_aug_state(hs, aug_state, u, times_is_active[k], tsteps[k], maps_sys[k])

        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end

    return aug_state_traj, u_traj
end

"""Checks if the final goal is reached
"""
function reached(specs::ProblemSpecs, aug_state)
    (x, t, k) = aug_state
    # Checks if the current mode is part of the target modes
    idx = findfirst(==(k), specs.Ns_target)
    if isnothing(idx)
        return false
    end
    # Checks spatial and temporal inclusion
    X_target = specs.Xs_target[idx]
    T_target = specs.Ts_target[idx]
    in_X = x ∈ X_target
    in_T = t ≥ T_target.lb[1] && t ≤ T_target.ub[1]
    return in_X && in_T
end

end
