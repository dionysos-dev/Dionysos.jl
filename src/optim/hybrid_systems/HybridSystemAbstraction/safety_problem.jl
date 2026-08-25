"""
    OptimizerSafetyProblem{T} <: Dionysos.Optim.AbstractLiftedControlOptimizer

Safety sub-solver of the hybrid family: lifts a
[`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) over a hybrid system onto the
[`HybridSymbolicModel`](@ref Dionysos.Symbolic.HybridSymbolicModel) abstraction and computes the
maximal controlled-invariant set on it.

The safe set is mode-indexed ([`HybridSpec`](@ref Dionysos.Problem.HybridSpec)); a mode the
specification says nothing about is bounded by its own state set rather than forbidden.

Set `"concrete_problem"` and `"abstract_system"`; read back `"abstract_controller"`,
`"invariant_set"` and `"invariant_set_complement"`.
"""
mutable struct OptimizerSafetyProblem{T} <: AbstractLiftedControlOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.HybridSymbolicModel}
    print_level::Int

    # Outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerSafetyProblem}
    abstract_problem::Union{Nothing, PR.SafetyProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    invariant_set::Any
    invariant_set_complement::Any
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

function reset!(model::OptimizerSafetyProblem)
    model.abstract_optimizer = nothing
    model.abstract_problem = nothing
    model.abstract_controller = nothing
    model.invariant_set = nothing
    model.invariant_set_complement = nothing
    model.abstract_problem_time_sec = 0.0
    model.success = false
    return model
end

abstract_optimizer_type(::OptimizerSafetyProblem) = OPDS.OptimizerSafetyProblem

function configure_abstract_optimizer!(model::OptimizerSafetyProblem, ::Any)
    model.print_level >= 1 && println("compute_largest_invariant_set! started")
    return
end

function extract_results!(model::OptimizerSafetyProblem, abstract_optimizer)
    model.invariant_set = abstract_optimizer.invariant_set
    model.invariant_set_complement = abstract_optimizer.invariant_set_complement
    model.print_level >= 1 &&
        model.invariant_set !== nothing &&
        println("Invariant set size: $(length(model.invariant_set))")
    return
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.HybridSymbolicModel,
)
    abstract_initial_set =
        _abstract_initial_states(abstract_system, concrete_problem.initial_set)
    _check_initial_nonempty(abstract_initial_set)

    abstract_safe_set = SY.states_satisfying(abstract_system, concrete_problem.safe_set)
    _check_nonempty(abstract_safe_set, "safe")

    return PR.SafetyProblem(
        SY.get_automaton(abstract_system),
        abstract_initial_set,
        abstract_safe_set,
        concrete_problem.time,
    )
end

function safe(concrete_problem::PR.SafetyProblem, aug_state)
    return PR.satisfies(concrete_problem.safe_set, aug_state...)
end
