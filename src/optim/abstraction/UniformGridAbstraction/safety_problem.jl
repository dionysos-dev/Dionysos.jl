"""
    OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer

An optimizer for solving **safety control problems** over symbolic system abstractions.

This solver takes as input a [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) and a symbolic abstraction of the system (e.g., a [`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList)), and computes a controller that ensures the system remains within a safe set over a time horizon or indefinitely.

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
  An instance of [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) that specifies the system, initial set, safe set, and horizon.

- `abstract_system` (**required**):  
  A symbolic abstraction of the system, e.g., obtained from [`OptimizerEmptyProblem`](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerEmptyProblem).

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
using Dionysos, JuMP
optimizer = MOI.instantiate(Dionysos.Optim.OptimizerSafetyProblem.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), my_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer)

time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
```
"""
mutable struct OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, Dionysos.Problem.SafetyProblem}
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}

    # Constructed parameters
    abstract_problem::Union{Nothing, Dionysos.Problem.SafetyProblem}
    abstract_controller::Union{Nothing, Dionysos.System.SymbolicController}
    abstract_problem_time_sec::T

    # Problem/Solver-Specific parameters
    invariant_set::Union{Nothing, Dionysos.Domain.DomainList}
    invariant_set_complement::Union{Nothing, Dionysos.Domain.DomainList}

    success::Bool
    print_level::Int
    function OptimizerSafetyProblem{T}() where {T}
        return new{T}(nothing, nothing, nothing, nothing, 0.0, nothing, nothing, false, 1)
    end
end

OptimizerSafetyProblem() = OptimizerSafetyProblem{Float64}()

MOI.is_empty(optimizer::OptimizerSafetyProblem) = optimizer.concrete_problem === nothing

function MOI.set(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerSafetyProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function MOI.optimize!(optimizer::OptimizerSafetyProblem)
    t_ref = time()

    # Ensure abstract system is set before proceeding
    if optimizer.abstract_system === nothing
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    end

    # Build the abstract problem from the concrete problem
    optimizer.abstract_problem =
        build_abstract_problem(optimizer.concrete_problem, optimizer.abstract_system)

    optimizer.print_level >= 1 && println("compute_controller_safe! started")
    # Compute the largest invariant set
    abstract_controller, invariant_set_symbols, invariant_set_complement_symbols =
        SY.compute_largest_invariant_set(
            optimizer.abstract_problem.system.autom,
            optimizer.abstract_problem.safe_set,
        )

    invariant_set = Dionysos.Symbolic.get_domain_from_states(
        optimizer.abstract_system,
        invariant_set_symbols,
    )
    invariant_set_complement = Dionysos.Symbolic.get_domain_from_states(
        optimizer.abstract_system,
        invariant_set_complement_symbols,
    )

    optimizer.abstract_controller = abstract_controller
    optimizer.invariant_set = invariant_set
    optimizer.invariant_set_complement = invariant_set_complement

    # Display results
    if ⊆(optimizer.abstract_problem.initial_set, invariant_set_symbols)
        optimizer.success = true
    end

    optimizer.print_level >= 1 && println("\n Safety: terminated with $(optimizer.success)")
    optimizer.abstract_problem_time_sec = time() - t_ref
    return
end

function build_abstract_problem(
    concrete_problem::Dionysos.Problem.SafetyProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
)
    return Dionysos.Problem.SafetyProblem(
        abstract_system,
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.initial_set,
            Dionysos.Domain.OUTER,
        ),
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.safe_set,
            Dionysos.Domain.INNER,
        ),
        concrete_problem.time, # TODO: This is continuous time, not the number of transitions
    )
end
