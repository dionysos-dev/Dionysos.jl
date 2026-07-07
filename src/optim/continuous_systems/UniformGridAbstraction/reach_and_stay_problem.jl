"""
    OptimizerReachAndStayProblem{T} <: MOI.AbstractOptimizer

An optimizer for solving **reach-and-stay control problems** over symbolic system abstractions.

This solver takes as input a [`ReachAndStayProblem`](@ref PR.ReachAndStayProblem) and a symbolic abstraction of the system (e.g., a [`SymbolicModel`](@ref SY.SymbolicModel)), and computes a controller that drives the system into the target set and keeps it there thereafter.

---

### Key Behavior

- Lifts the concrete reach-and-stay problem to the abstract domain and builds an `abstract_problem`.
- Computes the **winning set**, i.e., the largest set of abstract states from which the reach-and-stay specification can be enforced.
- Synthesizes an `abstract_controller` that guarantees the specification under worst-case transitions.
- The optimization is successful if `success == true` after calling `MOI.optimize!`.

---

### Parameters

#### Mandatory fields set by the user

- `concrete_problem` (**required**):  
  An instance of [`ReachAndStayProblem`](@ref PR.ReachAndStayProblem) specifying the system, initial set, target set, safe set, and time horizon (finite or infinite).

- `abstract_system` (**required**):  
  A symbolic abstraction of the system, typically produced by
  [`OptimizerAlternatingSimulationProblem`](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerAlternatingSimulationProblem).

#### Optional user-tunable parameters

- `early_stop` (optional, default = `false`):  
  Stops the fixed-point computation as soon as all abstract initial states belong to the winning set.
  This may significantly reduce computation time when only the feasibility of the given initial set is of interest.

- `print_level` (optional, default = `1`):  
  Controls verbosity:
    - `0`: silent
    - `1`: default (progress information)
    - `2`: verbose debug output

#### Internally computed fields

These fields are automatically filled in by `MOI.optimize!`.

- `abstract_problem`: The lifted reach-and-stay problem over the symbolic abstraction.
- `abstract_optimizer`: The discrete reach-and-stay optimizer used to solve the abstract problem.
- `abstract_controller`: A symbolic controller satisfying the reach-and-stay specification.
- `winning_set`: The maximal set of abstract states from which the specification can be enforced.
- `winning_set_complement`: The complement of the winning set.
- `success`: Boolean flag indicating whether all abstract initial states belong to the winning set.
- `abstract_problem_time_sec`: Time taken to solve the abstract reach-and-stay problem.
"""
mutable struct OptimizerReachAndStayProblem{T} <: MOI.AbstractOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_system::Any
    early_stop::Bool
    print_level::Int

    # outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerReachAndStayProblem}
    abstract_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    winning_set::Union{Nothing, MP.AbstractStateSet}
    winning_set_complement::Union{Nothing, MP.AbstractStateSet}
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerReachAndStayProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            false,
            1,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            zero(T),
        )
    end
end

OptimizerReachAndStayProblem() = OptimizerReachAndStayProblem{Float64}()

MOI.is_empty(optimizer::OptimizerReachAndStayProblem) =
    optimizer.concrete_problem === nothing

function MOI.set(
    model::OptimizerReachAndStayProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerReachAndStayProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerReachAndStayProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstract_problem(
    concrete_problem::PR.ReachAndStayProblem,
    abstract_system::SY.SymbolicModel,
)
    return PR.ReachAndStayProblem(
        SY.get_automaton(abstract_system),
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER),
        SY.get_states_from_set(abstract_system, concrete_problem.target_set, MP.INNER),
        SY.get_states_from_set(abstract_system, concrete_problem.safe_set, MP.INNER),
        concrete_problem.time,
    )
end

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing && error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    abstract_optimizer = MOI.instantiate(OPDS.OptimizerReachAndStayProblem)

    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)

    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("early_stop"),
        optimizer.early_stop,
    )

    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(abstract_optimizer)

    optimizer.abstract_optimizer = abstract_optimizer
    optimizer.abstract_controller = abstract_optimizer.controller

    optimizer.winning_set = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.winning_set),
    )

    optimizer.winning_set_complement = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.winning_set_complement),
    )

    optimizer.success = abstract_optimizer.success
    optimizer.abstract_problem_time_sec = time() - t0

    return
end
