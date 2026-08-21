mutable struct OptimizerSafetyProblem{T} <: OP.AbstractDionysosOptimizer
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

MOI.is_empty(optimizer::OptimizerSafetyProblem) = optimizer.concrete_problem === nothing

function MOI.get(model::OptimizerSafetyProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

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

function MOI.optimize!(optimizer::OptimizerSafetyProblem)
    t0 = time()

    optimizer.abstract_system === nothing && error("abstract_system not set")
    optimizer.concrete_problem === nothing && error("concrete_problem not set")

    abstract_system = optimizer.abstract_system

    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    optimizer.print_level >= 1 && println("compute_largest_invariant_set! started")

    abstract_optimizer = MOI.instantiate(OPDS.OptimizerSafetyProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(abstract_optimizer)

    optimizer.abstract_optimizer = abstract_optimizer
    optimizer.abstract_controller = abstract_optimizer.controller
    optimizer.invariant_set = abstract_optimizer.invariant_set
    optimizer.invariant_set_complement = abstract_optimizer.invariant_set_complement
    optimizer.success = abstract_optimizer.success

    if !isnothing(optimizer.invariant_set)
        optimizer.print_level >= 1 &&
            println("Invariant set size: $(length(optimizer.invariant_set))")
    end

    optimizer.abstract_problem_time_sec = time() - t0
    return
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.HybridSymbolicModel,
)
    concrete_initial_state = concrete_problem.initial_set
    q0 = SY.get_abstract_state(abstract_system, concrete_initial_state)
    q0 <= 0 && error("Initial augmented state is outside the abstract system domain.")
    abstract_initial_set = [q0]

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
