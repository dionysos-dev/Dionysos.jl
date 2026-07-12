"""
    OptimizerSafetyProblem{T} <: AbstractLiftedControlOptimizer

An optimizer for solving **safety control problems** over symbolic system abstractions.

This solver takes as input a [`SafetyProblem`](@ref PR.SafetyProblem) and a symbolic abstraction of the system (e.g., a [`SymbolicModelList`](@ref SY.SymbolicModelList)), and computes a controller that ensures the system remains within a safe set over a time horizon or indefinitely.

---

### Key Behavior

- Lifts the concrete safety problem to the abstract domain and builds an `abstract_problem`.
- Computes the **invariant set**, i.e., the largest set of abstract states from which all trajectories can be safely controlled.
- Synthesizes an `abstract_controller` that guarantees safety under worst-case transitions.
- The optimization is successful if `success == true` after calling `MOI.optimize!`.

---

### Parameters

#### Mandatory fields set by the user

- `concrete_problem` (**required**):  
  An instance of [`SafetyProblem`](@ref PR.SafetyProblem) that specifies the system, initial set, safe set, and horizon.

- `abstract_system` (**required**):  
  A symbolic abstraction of the system, e.g., obtained from [`OptimizerAlternatingSimulationProblem`](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerAlternatingSimulationProblem).

#### Optional user-tunable parameters

- `print_level` (optional, default = `1`):  
  Controls verbosity:
    - `0`: silent
    - `1`: default (info)
    - `2`: verbose debug output

#### Internally computed fields

These fields are automatically filled in by `MOI.optimize!`.

- `abstract_problem`: The lifted version of the safety problem in the symbolic domain.
- `abstract_problem_time_sec`: Time taken to solve the safety problem over the abstract system.
- `abstract_controller`: A controller mapping abstract states to admissible inputs that keep the system safe.
- `invariant_set`: The largest subset of abstract states from which safety can be maintained.
- `invariant_set_complement`: States from which safety cannot be guaranteed.
- `success`: Boolean flag indicating whether a valid invariant set and controller were found.

### Example

```julia
using Dionysos,
optimizer = MOI.instantiate(Dionysos.Optim.OptimizerSafetyProblem.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), my_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer)

time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
```
"""
mutable struct OptimizerSafetyProblem{T} <: AbstractLiftedControlOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.SafetyProblem}
    abstract_system::Any
    print_level::Int

    # outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerSafetyProblem}
    abstract_problem::Union{Nothing, PR.SafetyProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    invariant_set::Union{Nothing, MP.AbstractStateSet}
    invariant_set_complement::Union{Nothing, MP.AbstractStateSet}
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerSafetyProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
            1,       # print_level
            nothing, # abstract_optimizer
            nothing, # abstract_problem
            nothing, # abstract_controller
            nothing, # invariant_set
            nothing, # invariant_set_complement
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerSafetyProblem() = OptimizerSafetyProblem{Float64}()

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.SymbolicModel,
)
    return PR.SafetyProblem(
        SY.get_automaton(abstract_system),
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER),
        SY.get_states_from_set(abstract_system, concrete_problem.safe_set, MP.INNER),
        concrete_problem.time,
    )
end

abstract_optimizer_type(::OptimizerSafetyProblem) = OPDS.OptimizerSafetyProblem

function extract_results!(model::OptimizerSafetyProblem, abstract_optimizer)
    abs_sys = model.abstract_system
    model.invariant_set =
        SY.get_state_set_from_states(abs_sys, collect(abstract_optimizer.invariant_set))
    model.invariant_set_complement = SY.get_state_set_from_states(
        abs_sys,
        collect(abstract_optimizer.invariant_set_complement),
    )
    return
end
