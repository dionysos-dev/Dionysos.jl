mutable struct OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.TimedHybridSymbolicModel}
    print_level::Int

    # Outputs
    abstract_optimizer::Union{Nothing, SY.OptimizerSafetyProblem}
    abstract_problem::Union{Nothing, PR.SafetyProblem}
    abstract_controller::Union{Nothing, MS.AbstractMap}
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

function MOI.set(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerSafetyProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
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

    optimizer.abstract_system === nothing &&
        error("abstract_system not set")
    optimizer.concrete_problem === nothing &&
        error("concrete_problem not set")

    abstract_system = optimizer.abstract_system

    abstract_problem = build_abstract_problem(
        optimizer.concrete_problem,
        abstract_system,
    )
    optimizer.abstract_problem = abstract_problem

    optimizer.print_level >= 1 &&
        println("compute_largest_invariant_set! started")

    abstract_optimizer = MOI.instantiate(SY.OptimizerSafetyProblem)
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
    abstract_system::SY.TimedHybridSymbolicModel,
)
    concrete_initial_state = concrete_problem.initial_set
    abstract_initial_set = [SY.get_abstract_state(abstract_system, concrete_initial_state)]

    concrete_safe_set = concrete_problem.safe_set
    abstract_safe_set = SY.get_states_from_set(
        abstract_system,
        concrete_safe_set[1], # state sets
        concrete_safe_set[2], # time sets
        concrete_safe_set[3]; # mode list
        domain = MP.OUTER, # TODO MP.INNER
    )

    return PR.SafetyProblem(
        SY.get_automaton(abstract_system),
        abstract_initial_set,
        abstract_safe_set,
        concrete_problem.time,
    )
end

function safe(concrete_problem::PR.SafetyProblem, aug_state)
    (x, t, k) = aug_state
    (Xs_safe, Ts_safe, Ns_safe) = concrete_problem.safe_set

    idx = findfirst(==(k), Ns_safe)
    isnothing(idx) && return false

    X_set = Xs_safe[idx]
    T_set = Ts_safe[idx]

    in_X = x ∈ X_set
    in_T = T_set.lb[1] <= t <= T_set.ub[1]

    return in_X && in_T
end
