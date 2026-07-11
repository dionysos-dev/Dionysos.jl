"""
    OptimizerCoSafeLTLProblem{T} <: Dionysos.Optim.AbstractDionysosOptimizer

Abstraction-based solver for co-safe LTL control problems.

This optimizer:
1. lifts a concrete `CoSafeLTLProblem` to an abstract automaton problem,
2. calls the generic automaton-level co-safe LTL optimizer in `SY`,
3. stores the resulting abstract controller and solve status.
"""
mutable struct OptimizerCoSafeLTLProblem{T} <: AbstractLiftedControlOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_system::Any
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    abstract_optimizer::Union{Nothing, OPDS.OptimizerCoSafeLTLProblem}
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerCoSafeLTLProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
            false,   # sparse_input
            1,       # print_level
            nothing, # abstract_optimizer
            nothing, # abstract_problem
            nothing, # abstract_controller
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerCoSafeLTLProblem() = OptimizerCoSafeLTLProblem{Float64}()

function build_abstract_problem(
    concrete_problem::PR.CoSafeLTLProblem,
    abstract_system::SY.SymbolicModel,
)
    init_states =
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER)

    lab_abs = Dict{Symbol, Vector{Int}}()
    for (ap, setX) in concrete_problem.labeling
        sem = get(concrete_problem.ap_semantics, ap, MP.INNER)
        lab_abs[ap] = SY.get_states_from_set(abstract_system, setX, sem)
    end

    return PR.CoSafeLTLProblem(
        SY.get_automaton(abstract_system),
        init_states,
        concrete_problem.spec,
        lab_abs,
        concrete_problem.ap_semantics,
    )
end

abstract_optimizer_type(::OptimizerCoSafeLTLProblem) = OPDS.OptimizerCoSafeLTLProblem

function configure_abstract_optimizer!(model::OptimizerCoSafeLTLProblem, abstract_optimizer)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        model.sparse_input,
    )
    return
end
