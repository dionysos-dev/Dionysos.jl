"""
    OptimizerCoSafeLTLProblem{T} <: MOI.AbstractOptimizer

Abstraction-based solver for co-safe LTL control problems.

This optimizer:
1. lifts a concrete `CoSafeLTLProblem` to an abstract automaton problem,
2. calls the generic automaton-level co-safe LTL optimizer in `SY`,
3. stores the resulting abstract controller and solve status.
"""
mutable struct OptimizerCoSafeLTLProblem{T} <: MOI.AbstractOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    abstract_optimizer::Union{Nothing, SY.OptimizerCoSafeLTLProblem}
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_controller::Union{Nothing, MS.SystemWithOutput}
    qa0::Union{Nothing, Int}
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
            nothing, # qa0
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerCoSafeLTLProblem() = OptimizerCoSafeLTLProblem{Float64}()

MOI.is_empty(opt::OptimizerCoSafeLTLProblem) = opt.concrete_problem === nothing

function MOI.set(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute, v)
    return setproperty!(opt, Symbol(p.name), v)
end

function MOI.get(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute)
    return getproperty(opt, Symbol(p.name))
end

MOI.get(opt::OptimizerCoSafeLTLProblem, ::MOI.SolveTimeSec) = opt.abstract_problem_time_sec

# ============================================================
# Lifting: concrete CoSafeLTLProblem -> abstract CoSafeLTLProblem
# ============================================================

function build_abstract_problem(
    concrete_problem::PR.CoSafeLTLProblem,
    abstract_system::SY.SymbolicModelList,
)
    init_states =
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER)

    lab_abs = Dict{Symbol, Vector{Int}}()
    for (ap, setX) in concrete_problem.labeling
        sem = get(concrete_problem.ap_semantics, ap, MP.INNER)
        lab_abs[ap] = SY.get_states_from_set(abstract_system, setX, sem)
    end

    return PR.CoSafeLTLProblem(
        abstract_system.autom,
        init_states,
        concrete_problem.spec,
        lab_abs,
        concrete_problem.ap_semantics,
        concrete_problem.strict_spot,
    )
end

# ============================================================
# optimize!
# ============================================================

function MOI.optimize!(optimizer::OptimizerCoSafeLTLProblem)
    t0 = time()

    optimizer.abstract_system === nothing && error("abstract_system not set")
    optimizer.concrete_problem === nothing && error("concrete_problem not set")

    concrete_problem = optimizer.concrete_problem
    abs_sys = optimizer.abstract_system

    abstract_problem = build_abstract_problem(concrete_problem, abs_sys)
    optimizer.abstract_problem = abstract_problem

    abstract_optimizer = MOI.instantiate(SY.OptimizerCoSafeLTLProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        optimizer.sparse_input,
    )
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(abstract_optimizer)

    optimizer.abstract_optimizer = abstract_optimizer
    optimizer.abstract_controller = abstract_optimizer.controller
    optimizer.qa0 = abstract_optimizer.qa0
    optimizer.success = abstract_optimizer.success
    optimizer.abstract_problem_time_sec = time() - t0

    return
end