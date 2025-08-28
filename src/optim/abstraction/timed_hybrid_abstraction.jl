# """
# Module for timed hybrid system abstraction and controller synthesis.

# Provides tools for:
# - Optimal control problems on timed hybrid systems
# - Safety problems with invariant set computation  
# - Concrete controller synthesis from abstract controllers
# - Closed-loop trajectory simulation
# """
module TimedHybridAbstraction

import HybridSystems
import HybridSystems: HybridSystem
import MathematicalSystems
import StaticArrays: SVector
import Dionysos

# ================================================================
# Problem specification structures
# ================================================================

"""
    TimedHybridProblemSpecs{F}

Specification for timed hybrid control problems (optimal control or safety).

# Fields
- `initial_state::Tuple{AbstractVector{Float64}, Float64, Int}`: Initial augmented state ([x], t, mode_id)
- `Xs_target::Vector{<:Dionysos.Utils.HyperRectangle}`: Target/safe sets (spatial component)
- `Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}`: Target/safe sets (temporal component)  
- `Ns_target::Vector{Int}`: Target/safe mode indices
- `concret_cost_fun::F`: Concrete cost function for optimal control
- `problem_type::Symbol`: `:optimal_control` or `:safety`
- `time_horizon::Union{Real, Dionysos.Problem.Infinity}`: Time horizon constraint
"""
struct TimedHybridProblemSpecs{F} # Need to add initial region for safety problems (?)
    initial_state::Tuple{AbstractVector{Float64}, Float64, Int} # aug_state initial ([x], t, mode_id)
    Xs_target::Vector{<:Dionysos.Utils.HyperRectangle} # target sets or safe sets depending on problem type
    Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}
    Ns_target::Vector{Int}
    concret_cost_fun::F
    problem_type::Symbol # :optimal_control or :safety
    time_horizon::Union{Float64, Dionysos.Problem.Infinity} # User-specified time horizon
end

"""
    TimedHybridOptimalControlProblem(initial_state, Xs_target, Ts_target, Ns_target, cost_function, time_horizon)

Constructor for optimal control problems on timed hybrid systems.

# Arguments
- `initial_state`: Initial augmented state ([x], t, mode_id)
- `Xs_target`: Vector of spatial target sets (one per target mode)
- `Ts_target`: Vector of temporal target sets (one per target mode)  
- `Ns_target`: Vector of target mode indices
- `cost_function`: Cost function for transitions (aug_state, input) → cost
- `time_horizon`: Maximum time horizon (default: infinite)

# Returns
- `TimedHybridProblemSpecs`: Problem specification for optimal control
"""
function TimedHybridOptimalControlProblem(
    initial_state,
    Xs_target,
    Ts_target,
    Ns_target,
    cost_function,
    time_horizon = Dionysos.Problem.Infinity(),
)
    return TimedHybridProblemSpecs(
        initial_state,
        Xs_target,
        Ts_target,
        Ns_target,
        cost_function,
        :optimal_control,
        time_horizon,
    )
end

"""
    TimedHybridSafetyProblem(initial_state, Xs_safe, Ts_safe, Ns_safe, time_horizon)

Constructor for safety problems on timed hybrid systems.

# Arguments  
- `initial_state`: Initial augmented state ([x], t, mode_id)
- `Xs_safe`: Vector of spatial safe sets (one per safe mode)
- `Ts_safe`: Vector of temporal safe sets (one per safe mode)
- `Ns_safe`: Vector of safe mode indices  
- `time_horizon`: Maximum time horizon (default: infinite)

# Returns
- `TimedHybridProblemSpecs`: Problem specification for safety
"""
function TimedHybridSafetyProblem(
    initial_state,
    Xs_safe,
    Ts_safe,
    Ns_safe,
    time_horizon = Dionysos.Problem.Infinity(),
)
    # For safety problems, we use a dummy cost function 
    dummy_cost = (aug_state, u) -> 1.0
    return TimedHybridProblemSpecs(
        initial_state,
        Xs_safe,
        Ts_safe,
        Ns_safe,
        dummy_cost,
        :safety,
        time_horizon,
    )
end

# ================================================================
# Problem construction functions  
# ================================================================

# """
#     build_concrete_problem(problem_specs::TimedHybridProblemSpecs)

# Build concrete problem instance from problem specifications.

# # Arguments
# - `problem_specs`: Timed hybrid problem specifications

# # Returns  
# - `Dionysos.Problem.OptimalControlProblem` or `Dionysos.Problem.SafetyProblem`
# """
function build_concrete_problem(problem_specs::TimedHybridProblemSpecs)
    concrete_initial_set = problem_specs.initial_state

    if problem_specs.problem_type == :optimal_control
        concrete_target_set =
            (problem_specs.Xs_target, problem_specs.Ts_target, problem_specs.Ns_target)
        concrete_cost_fun = problem_specs.concret_cost_fun

        return Dionysos.Problem.OptimalControlProblem(
            problem_specs,  # system will be set in abstract problem
            concrete_initial_set,
            concrete_target_set,
            nothing,  # state_cost
            concrete_cost_fun,  # transition_cost
            problem_specs.time_horizon,
        )
    elseif problem_specs.problem_type == :safety
        concrete_safe_set = (
            problem_specs.Xs_target, # For safety, Xs_target are the safe sets
            problem_specs.Ts_target,
            problem_specs.Ns_target,
        )

        return Dionysos.Problem.SafetyProblem(
            problem_specs,  # system will be set in abstract problem
            concrete_initial_set,
            concrete_safe_set,
            problem_specs.time_horizon,
        )
    else
        error("Unknown problem type: $(problem_specs.problem_type)")
    end
end

# """
#     build_abstract_problem(concrete_problem, symmodel)

# Build abstract problem from concrete problem and symbolic model.

# # Arguments
# - `concrete_problem`: Concrete optimal control or safety problem
# - `symmodel`: Timed hybrid symbolic model

# # Returns
# - Abstract version of the problem with discrete state and input spaces
# """
function build_abstract_problem(
    concrete_problem::Union{
        Dionysos.Problem.OptimalControlProblem,
        Dionysos.Problem.SafetyProblem,
    },
    symmodel::Dionysos.Symbolic.SymbolicTimedHybridSystems.TimedHybridSymbolicModel,
)
    aug_state = concrete_problem.initial_set
    abstract_initial_set = [
        Dionysos.Symbolic.SymbolicTimedHybridSystems.get_abstract_state(
            symmodel,
            aug_state,
        ),
    ]

    if isa(concrete_problem, Dionysos.Problem.OptimalControlProblem)
        target_set = concrete_problem.target_set
        abstract_target_set =
            Dionysos.Symbolic.SymbolicTimedHybridSystems.get_states_from_set(
                symmodel,
                target_set[1],
                target_set[2],
                target_set[3],
            )

        return Dionysos.Problem.OptimalControlProblem(
            symmodel,
            abstract_initial_set,
            abstract_target_set,
            concrete_problem.state_cost,
            build_abstract_transition_cost(symmodel, concrete_problem.transition_cost),
            concrete_problem.time,
        )
    else # SafetyProblem
        safe_set = concrete_problem.safe_set
        # Use OUTER approximation for safe set to be more permissive
        abstract_safe_set =
            Dionysos.Symbolic.SymbolicTimedHybridSystems.get_states_from_set(
                symmodel,
                safe_set[1],
                safe_set[2],
                safe_set[3];
                domain = Dionysos.Domain.OUTER,
            )

        return Dionysos.Problem.SafetyProblem(
            symmodel,
            abstract_initial_set,
            abstract_safe_set,
            concrete_problem.time,
        )
    end
end
# ================================================================  
# Specialized invariant set computation for timed hybrid systems
# ================================================================

# """
# Compute the largest invariant set for timed hybrid automata systems.
# Unlike the standard algorithm, this version treats terminal states (states without outgoing transitions)
# as potentially safe, considering them as natural endpoints of the system evolution rather than unsafe states.

# # Arguments
# - `autom`: The automaton (temporal hybrid symbolic model)
# - `safelist`: Collection of initially safe states
# - `ControllerConstructor`: Function to create the controller

# # Returns
# - `controller`: Abstract controller
# - `invariant_set`: Set of states that can be kept safe indefinitely
# - `invariant_set_complement`: Set of states that cannot be kept safe
# """
function compute_largest_invariant_set_timed_hybrid( # (?) Will be changed, need to be discussed
    autom,
    safelist;
    ControllerConstructor::Function = () -> Dionysos.System.SymbolicControllerList(),
)
    controller = ControllerConstructor()
    nstates = Dionysos.Symbolic.get_n_state(autom)
    nsymbols = Dionysos.Symbolic.get_n_input(autom)

    # Initialize pairs table
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    # Compute valid (state, input) pairs
    for target in Dionysos.Symbolic.enum_states(autom)
        for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
            pairstable[source, symbol] = true
        end
    end

    # Count available actions per state
    nsymbolslist = sum(pairstable; dims = 2)

    # Initialize safe set
    safeset = Set(safelist)

    # For timed hybrid systems, we DON'T automatically remove terminal states
    # Instead, we classify them based on whether they represent valid system endpoints
    terminal_states = Set{Int}()
    for source in safeset
        if nsymbolslist[source] == 0
            push!(terminal_states, source)
        end
    end

    # Initialize unsafe set (everything not initially safe)
    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)

    # Disable transitions from initially unsafe states
    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end

    # Modified backward reachability analysis
    # We only remove states that can lead to genuinely unsafe states,
    # not terminal states that are naturally safe endpoints
    nextunsafeset = Set{Int}()
    iteration = 0

    while true
        iteration += 1
        # println(
        #     "   - Iteration $iteration: safe=$(length(safeset)), unsafe=$(length(unsafeset))",
        # )

        # Only propagate unsafety from genuinely unsafe states
        # (not from terminal states that are safe endpoints)
        truly_unsafe = setdiff(unsafeset, terminal_states)

        for target in truly_unsafe
            for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1

                    # A state becomes unsafe only if it has NO way to stay safe
                    # (not just because it can reach a terminal state)
                    if nsymbolslist[source] == 0 && !(source in terminal_states)
                        push!(nextunsafeset, source)
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        # Update sets
        setdiff!(safeset, nextunsafeset)
        union!(unsafeset, nextunsafeset)
        nextunsafeset = Set{Int}()
    end

    # println("   - Final safe states: $(length(safeset))")
    # println("   - Final unsafe states: $(length(unsafeset))")
    # println(
    #     "   - Terminal states kept safe: $(length(intersect(safeset, terminal_states)))",
    # )

    # Populate controller with valid actions from safe states
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                Dionysos.System.add_control!(controller, source, symbol)
            end
        end
    end

    # The complement is states that were initially safe but became unsafe
    invariant_set_complement = setdiff(Set(safelist), safeset)

    return controller, safeset, invariant_set_complement
end

# ================================================================
# Abstract problem solving
# ================================================================

# """
#     solve_abstract_problem(abstract_problem)

# Solve the abstract optimal control or safety problem.

# # Arguments
# - `abstract_problem`: Abstract problem instance  

# # Returns
# - `(abstract_controller, controllable_set)`: Controller and controllable states
# """
function solve_abstract_problem(
    abstract_problem::Union{
        Dionysos.Problem.OptimalControlProblem,
        Dionysos.Problem.SafetyProblem,
    },
)
    if isa(abstract_problem, Dionysos.Problem.OptimalControlProblem)
        abstract_controller, controllable_set_symbols, _, value_per_node =
            Dionysos.Symbolic.compute_worst_case_cost_controller(
                abstract_problem.system.symbolic_automaton,
                abstract_problem.target_set;
                initial_set = abstract_problem.initial_set,
                sparse_input = false,
                cost_function = abstract_problem.transition_cost,
            )

        println("value init_state : ", value_per_node[abstract_problem.initial_set[1]])

        if ⊆(abstract_problem.initial_set, controllable_set_symbols)
            println("✅ Optimal control problem is solvable: initial set is controllable")
        else
            println("⚠️ Warning: initial set is only partially controllable")
        end
    else # SafetyProblem
        println("\nSafe state number : ", length(collect(abstract_problem.safe_set)))
        println(
            "unsafe state number : ",
            Dionysos.Symbolic.get_n_state(abstract_problem.system.symbolic_automaton) -
            length(collect(abstract_problem.safe_set)),
        )

        # Use specialized timed hybrid invariant set computation
        abstract_controller, invariant_set_symbols, _ =
            compute_largest_invariant_set_timed_hybrid(
                abstract_problem.system.symbolic_automaton,
                collect(abstract_problem.safe_set);
                ControllerConstructor = () -> Dionysos.System.SymbolicControllerList(),
            )

        println("Controllable set size: $(length(invariant_set_symbols))")

        if ⊆(abstract_problem.initial_set, invariant_set_symbols)
            println("✅ Safety problem is solvable: initial set is safe-controllable")
        else
            println("⚠️ Warning: initial set is only partially safe-controllable")
        end

        # For safety problems, the controllable set is the invariant set
        controllable_set_symbols = invariant_set_symbols
    end

    return abstract_controller, controllable_set_symbols
end

# ================================================================  
# Concrete controller synthesis
# ================================================================

# """
#     solve_concrete_problem(symmodel, abstract_controller)

# Generate concrete controller from abstract controller.

# # Arguments
# - `symmodel`: Timed hybrid symbolic model
# - `abstract_controller`: Abstract controller

# # Returns
# - `Dionysos.System.BlackBoxContinuousController`: Concrete controller function
# """
function solve_concrete_problem(
    symmodel::Dionysos.Symbolic.SymbolicTimedHybridSystems.TimedHybridSymbolicModel,
    abstract_controller,
)
    function concrete_controller(aug_state)
        (x, t, k) = aug_state
        abstract_aug_state =
            Dionysos.Symbolic.SymbolicTimedHybridSystems.get_abstract_state(
                symmodel,
                aug_state,
            )

        if Dionysos.System.is_defined(abstract_controller, abstract_aug_state) == false
            println("⚠️ No action available for state $aug_state")
            return nothing
        end

        abstract_input =
            Dionysos.System.get_control(abstract_controller, abstract_aug_state)

        if Dionysos.Symbolic.SymbolicTimedHybridSystems.is_switching_input(
            symmodel.input_mapping,
            abstract_input,
        )
            transition_id = symmodel.input_mapping.global_to_switching[abstract_input]
            label = symmodel.input_mapping.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            return label
        else
            return Dionysos.Symbolic.SymbolicTimedHybridSystems.get_concrete_input(
                symmodel,
                abstract_input,
                k,
            )
        end
    end

    isdef = (aug_state) -> true # for the moment we assume the controller is defined for all augmented states
    return Dionysos.System.BlackBoxContinuousController(concrete_controller, isdef)
end

# ================================================================
# Abstract cost function construction
# ================================================================

# """
#     build_abstract_transition_cost(symmodel, concrete_cost_fun)

# Build abstract transition cost function from concrete cost function.

# # Arguments
# - `symmodel`: Timed hybrid symbolic model
# - `concrete_cost_fun`: Concrete cost function (aug_state, input) → cost

# # Returns  
# - Abstract cost function (state_int, input_int) → cost
# """
function build_abstract_transition_cost(
    symmodel::Dionysos.Symbolic.SymbolicTimedHybridSystems.TimedHybridSymbolicModel,
    concrete_cost_fun,
)
    function abstract_transition_cost(state_int, input_int)
        (x, t, k) = Dionysos.Symbolic.SymbolicTimedHybridSystems.get_concrete_state(
            symmodel,
            state_int,
        )
        aug_state = (x, t, k)
        if Dionysos.Symbolic.SymbolicTimedHybridSystems.is_switching_input(
            symmodel.input_mapping,
            input_int,
        )
            transition_id = symmodel.input_mapping.global_to_switching[input_int]
            label = symmodel.input_mapping.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            u = label
        else
            u = Dionysos.Symbolic.SymbolicTimedHybridSystems.get_concrete_input(
                symmodel,
                input_int,
                k,
            )
        end
        return concrete_cost_fun(aug_state, u)
    end
    return abstract_transition_cost
end

# ================================================================
# Main solving function
# ================================================================

"""
    solve_timed_hybrid_problem(hs, optimizer_factory_list, optimizer_kwargs_dict, problem_specs)

Solve a complete timed hybrid control problem (optimal control or safety).

# Arguments
- `hs::HybridSystem`: The hybrid system
- `optimizer_factory_list`: List of optimizer factories for each mode
- `optimizer_kwargs_dict`: Optimizer parameters for each mode  
- `problem_specs::TimedHybridProblemSpecs`: Problem specifications

# Returns
- Concrete controller for augmented states (x, t, mode)
"""
function solve_timed_hybrid_problem(
    hs::HybridSystem,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    problem_specs::TimedHybridProblemSpecs,
)
    hybrid_symmodel =
        Dionysos.Symbolic.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
            hs,
            optimizer_factory_list,
            optimizer_kwargs_dict,
        )

    concrete_problem = build_concrete_problem(problem_specs)

    abstract_problem = build_abstract_problem(concrete_problem, hybrid_symmodel)

    abstract_controller, _ = solve_abstract_problem(abstract_problem)

    concrete_controller = solve_concrete_problem(hybrid_symmodel, abstract_controller)

    return concrete_controller
end

# ================================================================
# Closed-loop simulation utilities  
# ================================================================

# """
#     get_next_aug_state(hs, aug_state, u, time_is_active, tstep, map_sys)

# Compute the next augmented state given current state and control input.

# # Arguments
# - `hs::HybridSystem`: The hybrid system
# - `aug_state`: Current augmented state (x, t, mode)
# - `u`: Control input (continuous or switching)
# - `time_is_active`: Whether time progresses in current mode
# - `tstep`: Time discretization step
# - `map_sys`: Discrete-time map for current mode

# # Returns
# - Next augmented state (x, t, mode)
# """
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
        # Explicit rounding of time to 10 decimals to avoid error propagation (?) ok or not 
        next_t = round(next_t; digits = 10)
        next_x = map_sys(x, u, tstep)
        return (next_x, next_t, k)
    end
end

"""
    get_closed_loop_trajectory(discretization_parameters, hs, problem_specs, controller, aug_state_0, nstep; stopping)

Generate a closed-loop trajectory using the synthesized controller.

# Arguments
- `discretization_parameters`: Time discretization parameters per mode  
- `hs::HybridSystem`: The hybrid system
- `problem_specs::TimedHybridProblemSpecs`: Problem specifications
- `controller`: Synthesized concrete controller
- `aug_state_0`: Initial augmented state
- `nstep`: Maximum number of simulation steps
- `stopping`: Stopping criterion function (default: never stop)

# Returns
- `(state_trajectory, control_trajectory)`: Trajectory and control sequence
"""
function get_closed_loop_trajectory(
    discretization_parameters::Vector{Tuple{Float64, Float64, Float64}},
    hs::HybridSystem,
    problem_specs::TimedHybridProblemSpecs,
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

# """
#     reached(specs, aug_state)

# Check if target is reached (optimal control) or safety is maintained (safety problems).

# # Arguments
# - `specs::TimedHybridProblemSpecs`: Problem specifications
# - `aug_state`: Current augmented state (x, t, mode)

# # Returns
# - `Bool`: True if goal is reached or safety is maintained
# """
function reached(specs::TimedHybridProblemSpecs, aug_state)
    (x, t, k) = aug_state
    # Checks if the current mode is part of the target/safe modes
    idx = findfirst(==(k), specs.Ns_target)
    if isnothing(idx)
        return false
    end
    # Checks spatial and temporal inclusion
    X_set = specs.Xs_target[idx]  # Target set for optimal control, safe set for safety
    T_set = specs.Ts_target[idx]
    in_X = x ∈ X_set
    in_T = t ≥ T_set.lb[1] && t ≤ T_set.ub[1]

    return in_X && in_T
end

end
