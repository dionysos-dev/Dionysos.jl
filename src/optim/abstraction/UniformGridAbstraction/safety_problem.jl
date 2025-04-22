"""
    OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer

An optimizer for solving **safety control problems** on symbolic system abstractions.

This solver takes a [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) as input along with a symbolic abstraction of the system (i.e., an [`abstract_system`](@ref Dionysos.Symbolic.SymbolicModelList)), and computes a **controller** that ensures all trajectories remain within a given safe set indefinitely (i.e., an **invariant set**).

### Key Behavior

- Lifts the concrete safety problem to the symbolic abstraction space (`abstract_system`) and constructs the corresponding `abstract_problem`.
- Computes the `invariant_set` — the largest set of abstract states from which safety can be guaranteed.
- Synthesizes a `abstract_controller` that keeps the abstract system within this invariant set under worst-case transitions.
- The solver is successful if the field `success` is `true` after `MOI.optimize!`.

### Parameters

#### Inputs

- `concrete_problem`: An instance of [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) that defines the safe set, system, and dynamics.
- `abstract_system`: The symbolic abstraction of the system, such as one produced by [`OptimizerEmptyProblem`](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerEmptyProblem).

#### Abstract Problem Fields

- `abstract_problem`: The lifted safety problem over the abstract model.
- `abstract_controller`: Abstract controller.
- `abstract_problem_time_sec`: Time to solve the abstract safety problem.

#### Result Sets

- `invariant_set`: The largest set of abstract states guaranteed to remain within the safe set.
- `invariant_set_complement`: The complement — states from which safety cannot be guaranteed under any control.

#### Miscellaneous

- `success`: A `Bool` flag that is `true` if a valid invariant set and controller were found.
- `print_level`: Controls verbosity:
    - `0`: silent
    - `1`: default (info)
    - `2`: verbose logging

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
    abstract_controller::Union{Nothing, Dionysos.Utils.SortedTupleSet{2, NTuple{2, Int}}}
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
        compute_largest_invariant_set(
            optimizer.abstract_problem.system,
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

function compute_largest_invariant_set(
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
    safelist,
)
    autom = abstract_system.autom
    contr = NewControllerList()
    nstates = autom.nstates
    nsymbols = autom.nsymbols
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable; dims = 2)

    # Remove unsafe states
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
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

    # Iterate until convergence
    while true
        for target in unsafeset
            for soursymb in Dionysos.Symbolic.pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end

    # Populate controller
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                Dionysos.Utils.push_new!(contr, (source, symbol))
            end
        end
    end
    unsafeset = setdiff(Set(safelist), safeset)
    return contr, safeset, unsafeset
end

function _compute_pairstable(pairstable, autom)
    for target in 1:(autom.nstates)
        for soursymb in Dionysos.Symbolic.pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end
