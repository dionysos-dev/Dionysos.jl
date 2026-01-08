module TimedHybridAbstraction
using MathOptInterface
using HybridSystems
import MathematicalSystems
MS = MathematicalSystems
using StaticArrays: SVector
using Dionysos

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

# ================================================================
# Problem specification structures
# ================================================================

# """
#     TimedHybridProblemSpecs{F}

# Specification for timed hybrid control problems (optimal control or safety).

# # Fields
# - `initial_state::Tuple{AbstractVector{Float64}, Float64, Int}`: Initial augmented state ([x], t, mode_id)
# - `Xs_target::Vector{<:Dionysos.Utils.HyperRectangle}`: Target/safe sets (spatial component)
# - `Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}`: Target/safe sets (temporal component)  
# - `Ns_target::Vector{Int}`: Target/safe mode indices
# - `concrete_cost_fun::F`: Concrete cost function for optimal control
# - `problem_type::Symbol`: `:optimal_control` or `:safety`
# - `time_horizon::Union{Real, Dionysos.Problem.Infinity}`: Time horizon constraint
# """
struct TimedHybridProblemSpecs{F}
    initial_state::Tuple{AbstractVector{Float64}, Float64, Int}
    Xs_target::Vector{<:Dionysos.Utils.HyperRectangle}
    Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}
    Ns_target::Vector{Int}
    concrete_cost_fun::F
    problem_type::Symbol
    time_horizon::Union{Float64, Dionysos.Problem.Infinity}
end

# """
#     TimedHybridOptimalControlProblem(initial_state, Xs_target, Ts_target, Ns_target, cost_function, time_horizon)

# Constructor for optimal control problems on timed hybrid systems.

# # Arguments
# - `initial_state`: Initial augmented state ([x], t, mode_id)
# - `Xs_target`: Vector of spatial target sets (one per target mode)
# - `Ts_target`: Vector of temporal target sets (one per target mode)  
# - `Ns_target`: Vector of target mode indices
# - `cost_function`: Cost function for transitions (aug_state, input) → cost
# - `time_horizon`: Maximum time horizon (default: infinite)

# # Returns
# - `TimedHybridProblemSpecs`: Problem specification for optimal control
# """
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

# """
#     TimedHybridSafetyProblem(initial_state, Xs_safe, Ts_safe, Ns_safe, time_horizon)

# Constructor for safety problems on timed hybrid systems.

# # Arguments  
# - `initial_state`: Initial augmented state ([x], t, mode_id)
# - `Xs_safe`: Vector of spatial safe sets (one per safe mode)
# - `Ts_safe`: Vector of temporal safe sets (one per safe mode)
# - `Ns_safe`: Vector of safe mode indices  
# - `time_horizon`: Maximum time horizon (default: infinite)

# # Returns
# - `TimedHybridProblemSpecs`: Problem specification for safety
# """
function TimedHybridSafetyProblem(
    initial_state,
    Xs_safe,
    Ts_safe,
    Ns_safe,
    time_horizon = Dionysos.Problem.Infinity(),
)
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
# Optimizer struct
# ================================================================

# """
#     Optimizer{T} <: MOI.AbstractOptimizer

# Abstraction-based solver for timed hybrid control problems using MathOptInterface.

# # Fields
# - `problem_specs`: Timed hybrid problem specifications
# - `concrete_problem`: Concrete optimal control or safety problem
# - `abstract_problem`: Abstract problem with discrete state/input spaces
# - `hybrid_system`: The hybrid system
# - `symbolic_model`: Timed hybrid symbolic model
# - `abstract_controller`: Abstract controller
# - `concrete_controller`: Concrete controller for augmented states
# - `optimizer_factory_list`: List of optimizer factories for each mode
# - `optimizer_kwargs_dict`: Optimizer parameters for each mode
# - `max_iterations`: Maximum number of iterations for solver
# - `solve_time_sec`: Time taken to solve the problem
# """
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    problem_specs::Union{Nothing, TimedHybridProblemSpecs}
    concrete_problem::Union{Nothing, Union{PR.OptimalControlProblem, PR.SafetyProblem}}
    abstract_problem::Union{Nothing, Union{PR.OptimalControlProblem, PR.SafetyProblem}}
    hybrid_system::Union{Nothing, HybridSystem}
    symbolic_model::Union{Nothing, SY.SymbolicTimedHybridSystems.TimedHybridSymbolicModel}
    abstract_controller::Union{Nothing, MS.AbstractMap}
    concrete_controller::Union{Nothing, MS.AbstractMap}
    optimizer_factory_list::Union{Nothing, Any}
    optimizer_kwargs_dict::Union{Nothing, Any}
    max_iterations::Union{Nothing, Int}
    solve_time_sec::T

    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0.0,
        )
    end
end

Optimizer() = Optimizer{Float64}()

# MOI Interface Methods
function MOI.is_empty(optimizer::Optimizer)
    return optimizer.problem_specs === nothing
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec)
    return optimizer.solve_time_sec
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "problem_specs"
        if !(value isa TimedHybridProblemSpecs)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    elseif param.name == "hybrid_system"
        if !(value isa HybridSystem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

# Set optimizer parameters
function set_optimizer!(
    optimizer::Optimizer,
    problem_specs::TimedHybridProblemSpecs,
    hybrid_system::HybridSystem,
    optimizer_factory_list,
    optimizer_kwargs_dict;
    max_iterations::Int = 1000,
)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("problem_specs"), problem_specs)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("hybrid_system"), hybrid_system)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_factory_list"),
        optimizer_factory_list,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        optimizer_kwargs_dict,
    )
    return MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iterations"), max_iterations)
end

# Optimize method
function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Ensure optimizer is properly configured
    if MOI.is_empty(optimizer)
        throw(
            MOI.InvalidStateException(
                "Optimizer is empty. Set problem_specs and hybrid_system.",
            ),
        )
    end

    # Extract components
    problem_specs = optimizer.problem_specs
    hs = optimizer.hybrid_system
    optimizer_factory_list = optimizer.optimizer_factory_list
    optimizer_kwargs_dict = optimizer.optimizer_kwargs_dict

    # Build symbolic model
    optimizer.symbolic_model =
        SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
            hs,
            optimizer_factory_list,
            optimizer_kwargs_dict,
        )

    # Build concrete and abstract problems
    optimizer.concrete_problem = build_concrete_problem(problem_specs)
    optimizer.abstract_problem =
        build_abstract_problem(optimizer.concrete_problem, optimizer.symbolic_model)

    # Solve abstract problem
    optimizer.abstract_controller, _ = solve_abstract_problem(optimizer.abstract_problem)

    # Synthesize concrete controller
    optimizer.concrete_controller =
        solve_concrete_problem(optimizer.symbolic_model, optimizer.abstract_controller)

    return optimizer.solve_time_sec = time() - t_ref
end

# ================================================================
# Problem construction functions
# ================================================================

function build_concrete_problem(problem_specs::TimedHybridProblemSpecs)
    concrete_initial_set = problem_specs.initial_state

    if problem_specs.problem_type == :optimal_control
        concrete_target_set =
            (problem_specs.Xs_target, problem_specs.Ts_target, problem_specs.Ns_target)
        concrete_cost_fun = problem_specs.concrete_cost_fun

        return Dionysos.Problem.OptimalControlProblem(
            problem_specs,
            concrete_initial_set,
            concrete_target_set,
            nothing,
            concrete_cost_fun,
            problem_specs.time_horizon,
        )
    elseif problem_specs.problem_type == :safety
        concrete_safe_set =
            (problem_specs.Xs_target, problem_specs.Ts_target, problem_specs.Ns_target)

        return Dionysos.Problem.SafetyProblem(
            problem_specs,
            concrete_initial_set,
            concrete_safe_set,
            problem_specs.time_horizon,
        )
    else
        error("Unknown problem type: $(problem_specs.problem_type)")
    end
end

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
    else
        safe_set = concrete_problem.safe_set
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
mutable struct SymbolicControlTable
    U::Vector{Vector{Int}}
end
SymbolicControlTable(nstates::Int) = SymbolicControlTable([Int[] for _ in 1:nstates])

function ensure_states!(C::SymbolicControlTable, nstates::Int)
    n = length(C.U)
    nstates <= n && return
    resize!(C.U, nstates)
    for i in (n + 1):nstates
        C.U[i] = Int[]
    end
end

add_control!(C::SymbolicControlTable, qs::Int, u::Int) =
    (ensure_states!(C, qs); push!(C.U[qs], u))
is_defined(C::SymbolicControlTable, qs::Int) = !isempty(C.U[qs])

struct PredicateDomain{F}
    ;
    pred::F;
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function to_ms_abstract_controller(C::SymbolicControlTable)
    qfun = (qs::Int) -> C.U[qs]
    X = PredicateDomain((qs::Int) -> is_defined(C, qs))
    return MS.ConstrainedBlackBoxMap(1, 1, qfun, X)  # (qs)->Vector{Int}
end

function compute_largest_invariant_set_timed_hybrid(autom, safelist)
    nstates = Dionysos.Symbolic.get_n_state(autom)
    nsymbols = Dionysos.Symbolic.get_n_input(autom)

    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    for target in Dionysos.Symbolic.enum_states(autom)
        for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
            pairstable[source, symbol] = true
        end
    end

    nsymbolslist = sum(pairstable; dims = 2)

    safeset = Set(safelist)
    terminal_states = Set{Int}()
    for source in safeset
        if nsymbolslist[source] == 0
            push!(terminal_states, source)
        end
    end

    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)

    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end

    nextunsafeset = Set{Int}()
    iteration = 0
    controller = SymbolicControlTable(nstates)

    while true
        iteration += 1
        truly_unsafe = setdiff(unsafeset, terminal_states)

        for target in truly_unsafe
            for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1

                    if nsymbolslist[source] == 0 && !(source in terminal_states)
                        push!(nextunsafeset, source)
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        union!(unsafeset, nextunsafeset)
        nextunsafeset = Set{Int}()
    end

    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                add_control!(controller, source, symbol)
            end
        end
    end

    invariant_set_complement = setdiff(Set(safelist), safeset)
    return to_ms_abstract_controller(controller), safeset, invariant_set_complement
end

# ================================================================
# Abstract problem solving
# ================================================================

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
    else
        println("\nSafe state number : ", length(collect(abstract_problem.safe_set)))
        println(
            "unsafe state number : ",
            Dionysos.Symbolic.get_n_state(abstract_problem.system.symbolic_automaton) -
            length(collect(abstract_problem.safe_set)),
        )

        abstract_controller, invariant_set_symbols, _ =
            compute_largest_invariant_set_timed_hybrid(
                abstract_problem.system.symbolic_automaton,
                collect(abstract_problem.safe_set),
            )

        println("Controllable set size: $(length(invariant_set_symbols))")

        if ⊆(abstract_problem.initial_set, invariant_set_symbols)
            println("✅ Safety problem is solvable: initial set is safe-controllable")
        else
            println("⚠️ Warning: initial set is only partially safe-controllable")
        end

        controllable_set_symbols = invariant_set_symbols
    end

    return abstract_controller, controllable_set_symbols
end

# ================================================================
# Concrete controller synthesis
# ================================================================

function solve_concrete_problem(
    symmodel::Dionysos.Symbolic.SymbolicTimedHybridSystems.TimedHybridSymbolicModel,
    abstract_controller::MS.AbstractMap,
)
    k_abs = abstract_controller.h  # q -> Vector{Int} or Int or nothing

    function f(aug_state)
        q = Dionysos.Symbolic.SymbolicTimedHybridSystems.get_abstract_state(
            symmodel,
            aug_state,
        )

        us = k_abs(q)
        us === nothing && return nothing
        (us isa AbstractVector && isempty(us)) && return nothing

        u_abs = us isa AbstractVector ? first(us) : us

        (_, _, k) = aug_state
        if Dionysos.Symbolic.SymbolicTimedHybridSystems.is_switching_input(
            symmodel.input_mapping,
            u_abs,
        )
            transition_id = symmodel.input_mapping.global_to_switching[u_abs]
            return symmodel.input_mapping.switch_labels[transition_id]  # string label
        else
            return Dionysos.Symbolic.SymbolicTimedHybridSystems.get_concrete_input(
                symmodel,
                u_abs,
                k,
            )
        end
    end

    X = PredicateDomain(_ -> true)  # or add a real domain check if you want
    nx = 1  # “dimension” is not super meaningful for tuples; keep consistent with MS expectations in your codebase
    nu = 1
    return MS.ConstrainedBlackBoxMap(nx, nu, f, X)
end

# ================================================================
# Abstract cost function construction
# ================================================================

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
            label = symmodel.input_mapping.switch_labels[transition_id]
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

# """
#     solve_timed_hybrid_problem(hs, optimizer_factory_list, optimizer_kwargs_dict, problem_specs; max_iterations)

# Solve a complete timed hybrid control problem using the MOI optimizer.

# # Arguments
# - `hs::HybridSystem`: The hybrid system
# - `optimizer_factory_list`: List of optimizer factories for each mode
# - `optimizer_kwargs_dict`: Optimizer parameters for each mode  
# - `problem_specs::TimedHybridProblemSpecs`: Problem specifications
# - `max_iterations`: Maximum number of iterations (default: 1000)

# # Returns
# - `Dionysos.System.BlackBoxContinuousController`: Concrete controller
# """
function solve_timed_hybrid_problem(
    hs::HybridSystem,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    problem_specs::TimedHybridProblemSpecs;
    max_iterations::Int = 1000,
)
    optimizer = Optimizer()
    set_optimizer!(
        optimizer,
        problem_specs,
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict;
        max_iterations = max_iterations,
    )
    MOI.optimize!(optimizer)
    return optimizer.concrete_controller
end

# ================================================================
# Closed-loop simulation utilities
# ================================================================

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

        reset_map = HybridSystems.resetmap(hs, transition)
        augmented_source_state = vcat(x, t)
        reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
        next_x = reset_result[1:(end - 1)]
        next_t = reset_result[end]
        next_k = target_mode
        next_t = round(next_t; digits = 10)
        return (next_x, next_t, next_k)
    else
        next_t = time_is_active ? t + tstep : 0.0
        next_t = round(next_t; digits = 10)
        next_x = map_sys(x, u, tstep)
        return (next_x, next_t, k)
    end
end

function get_closed_loop_trajectory(
    discretization_parameters::Vector{Tuple{Float64, Float64, Float64}},
    hs::HybridSystem,
    problem_specs::TimedHybridProblemSpecs,
    controller,
    aug_state_0,
    nstep;
    stopping = (x) -> false,
)
    kmap = controller.h
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
        u = kmap(aug_state)
        u === nothing && break

        (_, _, k) = aug_state
        aug_state =
            get_next_aug_state(hs, aug_state, u, times_is_active[k], tsteps[k], maps_sys[k])

        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end

    return aug_state_traj, u_traj
end

function reached(specs::TimedHybridProblemSpecs, aug_state)
    (x, t, k) = aug_state
    idx = findfirst(==(k), specs.Ns_target)
    if isnothing(idx)
        return false
    end
    X_set = specs.Xs_target[idx]
    T_set = specs.Ts_target[idx]
    in_X = x ∈ X_set
    in_T = t ≥ T_set.lb[1] && t ≤ T_set.ub[1]

    return in_X && in_T
end

end # module