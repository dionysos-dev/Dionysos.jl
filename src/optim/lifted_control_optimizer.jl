# ------------------------------------------------------------
# Shared template for every optimizer that lifts a concrete control problem onto
# a symbolic abstraction.
#
# Each one builds the abstract problem, runs the matching discrete-system
# optimizer, then reads results back. Only the discrete optimizer type, the
# options forwarded to it, and the extracted result sets differ per problem —
# those are the three hooks.
#
# It lives in `Optim` rather than in one solver family because both the
# uniform-grid family and the hybrid family lift problems this way. The hooks are
# therefore `Optim`-level generic functions, and a family adds methods to them.
# ------------------------------------------------------------

"""
    AbstractLiftedControlOptimizer <: AbstractDionysosOptimizer

Supertype for a solver that lifts a concrete control problem onto a symbolic abstraction and
delegates the synthesis to a discrete-system optimizer.

A subtype needs the fields `concrete_problem`, `abstract_system`, `abstract_problem`,
`abstract_optimizer`, `abstract_controller`, `success`, `print_level` and
`abstract_problem_time_sec`, and methods for:

- [`build_abstract_problem`](@ref)`(concrete_problem, abstract_system)` — lift the problem;
- [`abstract_optimizer_type`](@ref)`(model)` — the discrete optimizer to run;
- [`configure_abstract_optimizer!`](@ref)`(model, abstract_optimizer)` — forward extra options
  (optional);
- [`extract_results!`](@ref)`(model, abstract_optimizer)` — copy back the result sets (optional).

`MOI.is_empty`, `MOI.get(::MOI.SolveTimeSec)` and `MOI.optimize!` then come for free.
"""
abstract type AbstractLiftedControlOptimizer <: AbstractDionysosOptimizer end

MOI.is_empty(model::AbstractLiftedControlOptimizer) = model.concrete_problem === nothing

MOI.get(model::AbstractLiftedControlOptimizer, ::MOI.SolveTimeSec) =
    model.abstract_problem_time_sec

"""
    build_abstract_problem(concrete_problem, abstract_system) -> Problem.ProblemType

Lift a concrete control problem onto `abstract_system`. One method per (problem, abstraction)
pair, defined next to the optimizer that runs it.
"""
function build_abstract_problem end

"The discrete-system optimizer type an [`AbstractLiftedControlOptimizer`](@ref) delegates to."
abstract_optimizer_type(model::AbstractLiftedControlOptimizer) =
    error("abstract_optimizer_type not implemented for $(typeof(model))")

"Forward the options only this problem understands to the discrete optimizer. Optional."
configure_abstract_optimizer!(::AbstractLiftedControlOptimizer, ::Any) = nothing

"Copy the result sets only this problem produces back off the discrete optimizer. Optional."
extract_results!(::AbstractLiftedControlOptimizer, ::Any) = nothing

function MOI.optimize!(model::AbstractLiftedControlOptimizer)
    t0 = time()

    model.abstract_system === nothing && error("abstract_system not set")
    model.concrete_problem === nothing && error("concrete_problem not set")

    abstract_problem = build_abstract_problem(model.concrete_problem, model.abstract_system)
    model.abstract_problem = abstract_problem

    abstract_optimizer = MOI.instantiate(abstract_optimizer_type(model))
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("print_level"), model.print_level)
    configure_abstract_optimizer!(model, abstract_optimizer)

    MOI.optimize!(abstract_optimizer)

    model.abstract_optimizer = abstract_optimizer
    model.abstract_controller = abstract_optimizer.controller
    model.success = abstract_optimizer.success
    extract_results!(model, abstract_optimizer)

    model.abstract_problem_time_sec = time() - t0
    return
end
