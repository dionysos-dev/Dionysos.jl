mutable struct OptimizerOptimalControlProblem{T} <: OP.AbstractDionysosOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_transition_cost::Union{Nothing, Any}

    # Solver parameters
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # Outputs / internals
    abstract_optimizer::Union{Nothing, OPDS.OptimizerOptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    controllable_set::Union{Nothing, MP.AbstractStateSet}
    uncontrollable_set::Union{Nothing, MP.AbstractStateSet}
    value_fun_tab::Any
    abstract_value_function::Any
    concrete_value_function::Any
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
            nothing, # abstract_transition_cost
            false,   # early_stop
            false,   # sparse_input
            1,       # print_level
            nothing, # abstract_optimizer
            nothing, # abstract_problem
            nothing, # abstract_controller
            nothing, # controllable_set
            nothing, # uncontrollable_set
            nothing, # value_fun_tab
            nothing, # abstract_value_function
            nothing, # concrete_value_function
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerOptimalControlProblem() = OptimizerOptimalControlProblem{Float64}()

MOI.is_empty(optimizer::OptimizerOptimalControlProblem) =
    optimizer.concrete_problem === nothing

function MOI.get(model::OptimizerOptimalControlProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function reset!(model::OptimizerOptimalControlProblem)
    model.abstract_optimizer = nothing
    model.abstract_problem = nothing
    model.abstract_controller = nothing
    model.abstract_problem_time_sec = 0.0

    model.controllable_set = nothing
    model.uncontrollable_set = nothing
    model.value_fun_tab = nothing
    model.abstract_value_function = nothing
    model.concrete_value_function = nothing

    model.success = false
    return model
end

function build_abstract_value_function(value_fun_tab)
    return state -> value_fun_tab[state]
end

function build_concrete_value_function(
    abstract_system,
    abstract_value_function;
    default_value = Inf,
)
    function concrete_value_function(x)
        states = SY.get_abstract_states(abstract_system, x)
        isempty(states) && return default_value

        vals = Float64[]
        for q in states
            v = abstract_value_function(q)
            isfinite(v) && push!(vals, v)
        end

        isempty(vals) && return default_value
        return minimum(vals)
    end
    return concrete_value_function
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t0 = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing && error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    abstract_problem = build_abstract_problem(
        optimizer.concrete_problem,
        abstract_system,
        optimizer.abstract_transition_cost,
    )
    optimizer.abstract_problem = abstract_problem

    optimizer.print_level >= 1 && println("compute_controller_reachability! started")

    abstract_optimizer = MOI.instantiate(OPDS.OptimizerOptimalControlProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("early_stop"),
        optimizer.early_stop,
    )
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

    optimizer.controllable_set = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.controllable_set),
    )
    optimizer.uncontrollable_set = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.uncontrollable_set),
    )

    optimizer.value_fun_tab = abstract_optimizer.value_fun_tab
    optimizer.abstract_value_function =
        build_abstract_value_function(optimizer.value_fun_tab)
    optimizer.concrete_value_function =
        build_concrete_value_function(abstract_system, optimizer.abstract_value_function)

    optimizer.success = abstract_optimizer.success

    optimizer.print_level >= 1 &&
        println("\n Reachability: terminated with $(optimizer.success)")

    optimizer.abstract_problem_time_sec = time() - t0
    return
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.SymbolicModelList,
    abstract_transition_cost,
)
    @warn("The `state_cost` is not yet fully implemented")

    # `INNER` for the safe set: a cell only partly inside it contains unsafe points, so it
    # cannot be certified safe.
    concrete_safe_set = concrete_problem.safe_set
    abstract_safe_set =
        concrete_safe_set === nothing ? nothing :
        SY.get_states_from_set(abstract_system, concrete_safe_set, MP.INNER)

    return PR.OptimalControlProblem(
        SY.get_automaton(abstract_system),
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER),
        SY.get_states_from_set(abstract_system, concrete_problem.target_set, MP.INNER),
        concrete_problem.state_cost,
        abstract_transition_cost,
        concrete_problem.time;
        safe_set = abstract_safe_set,
    )
end
