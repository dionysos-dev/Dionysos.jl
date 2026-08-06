"""
    OptimizerOptimalControlProblem{T} <: AbstractLiftedControlOptimizer

An optimizer that solves reachability or reach-avoid **optimal control problems** using symbolic abstractions of the system.

This solver takes as input a concrete problem (typically an instance of [`OptimalControlProblem`](@ref PR.OptimalControlProblem)) and a symbolic abstraction of the system (i.e., an [`abstract_system`](@ref SY.SymbolicModel)). It then solves the **abstract** version of the control problem.

### Key Behavior

- Lifts the concrete problem to the symbolic abstraction space (`abstract_system`) and constructs the corresponding `abstract_problem`.
- Computes the `controllable_set` — the largest set of abstract states from which reachability can be guaranteed.
- Synthesizes an `abstract_controller` that brings the system to the target set under worst-case dynamics.
- Computes the `abstract_value_function` that maps each state (cell) to the worst-case number of steps needed to reach the target.
- The solver is successful if the field `success` is `true` after `MOI.optimize!`.

---

### Parameters

#### Mandatory fields set by the user

- `concrete_problem` (**required**):  
  An instance of [`OptimalControlProblem`](@ref PR.OptimalControlProblem) that defines the reach-avoid task (system, initial set, target, costs, horizon).

- `abstract_system` (**required**):  
  The symbolic abstraction of the system, usually obtained from an abstraction optimizer such as [`OptimizerAlternatingSimulationProblem`](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerAlternatingSimulationProblem).

#### Optional user-tunable parameters

- `early_stop` (optional, default = `false`):  
  If `true`, the fixpoint algorithm stops early when the initial set is fully contained in the controllable set.  
  If `false`, it computes the entire maximal controllable set.

- `sparse_input` (optional, default = `false`):  
  If `true`, uses a sparse representation of the transition table, reducing memory usage when the number of inputs is large but only few are admissible per state (e.g., in [`determinized abstractions`](@ref SY.determinize_symbolic_model), with `new_input = (input, target)`).

- `print_level` (optional, default = `1`):  
  Controls verbosity:  
    - `0`: silent  
    - `1`: default  
    - `2`: detailed logging

#### Internally computed fields

These fields are generated automatically during `MOI.optimize!`.

- `abstract_problem`: The lifted version of the concrete problem over the abstract system.
- `abstract_problem_time_sec`: Time taken to solve the abstract problem.
- `abstract_controller`: A controller mapping abstract states to control inputs.
- `controllable_set`: Set of abstract states from which the target is reachable.
- `uncontrollable_set`: Complementary states with no admissible reachability strategy.
- `value_fun_tab`: Tabular value function over abstract states (e.g., cost-to-go or step count).
- `abstract_value_function`: Functional form of the abstract value function.
- `concrete_value_function`: Functional form of the value function mapped back to the original system.
- `success`: Boolean flag indicating whether the solver completed successfully.

### Example

```julia
using Dionysos, JuMP
optimizer = MOI.instantiate(Dionysos.Optim.OptimizerOptimalControlProblem.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), my_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer)

time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
abstract_value_function = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
concrete_value_function = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
```
"""
mutable struct OptimizerOptimalControlProblem{T} <: AbstractLiftedControlOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Any
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerOptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}

    controllable_set::Union{Nothing, MP.AbstractStateSet}
    uncontrollable_set::Union{Nothing, MP.AbstractStateSet}
    value_fun_tab::Any
    abstract_value_function::Any
    concrete_value_function::Any

    success::Bool
    abstract_problem_time_sec::T

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
            false,   # early_stop
            false,   # sparse_input
            1,       # print_level
            nothing, # abstract_optimizer
            nothing, # abstract_problem
            nothing, # abstract_controller
            nothing, # controllable_set
            nothing, # uncontrollable_set
            nothing, # value_fun_tab
            nothing, # abstract_value_function
            nothing, # concrete_value_function
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerOptimalControlProblem() = OptimizerOptimalControlProblem{Float64}()

function reset!(model::OptimizerOptimalControlProblem)
    model.abstract_problem = nothing
    model.abstract_system = nothing
    model.abstract_controller = nothing
    model.abstract_problem_time_sec = 0.0

    model.controllable_set = nothing
    model.uncontrollable_set = nothing

    model.value_fun_tab = nothing
    model.abstract_value_function = nothing
    model.concrete_value_function = nothing

    model.success = false
    return model
end

function build_concrete_value_function(abstract_system, abstract_value_function)
    function concrete_value_function(x)
        state = SY.get_abstract_state(abstract_system, x)
        return abstract_value_function(state)
    end
    return concrete_value_function
end

function get_abstract_transition_cost(abstract_system, concrete_transition_cost)
    concrete_transition_cost === nothing && return nothing

    function abstract_transition_cost(q, s)
        x = SY.get_concrete_state(abstract_system, q)
        u = SY.get_concrete_input(abstract_system, s)
        return concrete_transition_cost(x, u)
    end

    return abstract_transition_cost
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.AbstractSymbolicModel,
)
    @warn("The `state_cost` is not yet fully implemented")

    return PR.OptimalControlProblem(
        SY.get_automaton(abstract_system),
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER),
        SY.get_states_from_set(abstract_system, concrete_problem.target_set, MP.INNER),
        concrete_problem.state_cost,
        get_abstract_transition_cost(abstract_system, concrete_problem.transition_cost),
        concrete_problem.time;
        safe_set = get_abstract_safe_set(abstract_system, concrete_problem.safe_set),
    )
end

"""
    get_abstract_safe_set(abstract_system, concrete_safe_set)

Lift the optional safe set of a reach-avoid problem to abstract states, keeping only cells
lying **entirely** inside it (`MP.INNER`) — a cell straddling the boundary contains unsafe
points, so it cannot be certified safe. `nothing` passes through, meaning the whole domain.
"""
function get_abstract_safe_set(abstract_system, concrete_safe_set)
    concrete_safe_set === nothing && return nothing
    return SY.get_states_from_set(abstract_system, concrete_safe_set, MP.INNER)
end

abstract_optimizer_type(::OptimizerOptimalControlProblem) =
    OPDS.OptimizerOptimalControlProblem

function configure_abstract_optimizer!(
    model::OptimizerOptimalControlProblem,
    abstract_optimizer,
)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("early_stop"), model.early_stop)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        model.sparse_input,
    )
    return
end

function extract_results!(model::OptimizerOptimalControlProblem, abstract_optimizer)
    abs_sys = model.abstract_system

    model.value_fun_tab = abstract_optimizer.value_fun_tab
    model.abstract_value_function = abstract_optimizer.value_function
    model.concrete_value_function =
        build_concrete_value_function(abs_sys, model.abstract_value_function)

    # The controllable/uncontrollable sets are grid-mapping-backed inspection sets;
    # a product model (e.g. clock-lifted) has no single mapping, so they stay `nothing`.
    if abs_sys isa SY.SymbolicModel
        model.controllable_set = SY.get_state_set_from_states(
            abs_sys,
            collect(abstract_optimizer.controllable_set),
        )
        model.uncontrollable_set = SY.get_state_set_from_states(
            abs_sys,
            collect(abstract_optimizer.uncontrollable_set),
        )
    end
    return
end
