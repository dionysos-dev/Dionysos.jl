# ------------------------------------------------------------
# Shared template for the four control-problem optimizers below.
#
# Each lifts a concrete problem to the abstract automaton, runs the matching
# discrete-system optimizer, then reads results back. Only the discrete
# optimizer type, the options forwarded to it, and the extracted result sets
# differ per problem — those are the three hooks.
# ------------------------------------------------------------

abstract type AbstractLiftedControlOptimizer <: OP.AbstractDionysosOptimizer end

MOI.is_empty(model::AbstractLiftedControlOptimizer) = model.concrete_problem === nothing

MOI.get(model::AbstractLiftedControlOptimizer, ::MOI.SolveTimeSec) =
    model.abstract_problem_time_sec

# Hooks:
#   abstract_optimizer_type(model)              -> the OPDS optimizer type to run
#   configure_abstract_optimizer!(model, ao)    -> forward extra options (default: none)
#   extract_results!(model, ao)                 -> copy back result sets (default: none)
# `build_abstract_problem(concrete, abstract_system)` is a per-problem method
# already defined next to each optimizer.
abstract_optimizer_type(model::AbstractLiftedControlOptimizer) =
    error("abstract_optimizer_type not implemented for $(typeof(model))")
configure_abstract_optimizer!(::AbstractLiftedControlOptimizer, abstract_optimizer) =
    nothing
extract_results!(::AbstractLiftedControlOptimizer, abstract_optimizer) = nothing

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
