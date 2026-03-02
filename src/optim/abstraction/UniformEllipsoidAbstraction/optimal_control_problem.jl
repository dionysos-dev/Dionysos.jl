mutable struct OptimizerOptimalControlProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem}
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}
    abstract_transition_cost::Union{Nothing, Any}

    # Common parameters
    abstract_problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem}
    abstract_controller::Union{Nothing, MS.ConstrainedBlackBoxMap}
    abstract_problem_time_sec::T

    # Specific parameters
    early_stop::Union{Nothing, Bool}
    sparse_input::Bool
    controllable_set::Union{Nothing, MP.AbstractStateSet}   
    uncontrollable_set::Union{Nothing, MP.AbstractStateSet} 
    value_fun_tab::Union{Nothing, Any} # Value function in tabular form, Inf means uncontrollable state
    abstract_value_function::Union{Nothing, Any}
    concrete_value_function::Union{Nothing, Any}

    success::Bool
    print_level::Int

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0.0,
            false,
            false,
            nothing,
            nothing,
            false,
            1,
        )
    end
end

OptimizerOptimalControlProblem() = OptimizerOptimalControlProblem{Float64}()

MOI.is_empty(optimizer::OptimizerOptimalControlProblem) =
    optimizer.concrete_problem === nothing

function MOI.set(
    model::OptimizerOptimalControlProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerOptimalControlProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerOptimalControlProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function reset!(model::OptimizerOptimalControlProblem)
    model.abstract_problem = nothing
    model.abstract_system = nothing
    model.abstract_transition_cost = nothing
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
    return abstract_value_function(state) = value_fun_tab[state]
end

function build_concrete_value_function(
    abstract_system,
    abstract_value_function;
    default_value = Inf
)
    function concrete_value_function(x)
        # get all abstract states covering x
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
    t_ref = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing &&
        error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    # Build abstract problem
    optimizer.abstract_problem =
        build_abstract_problem(optimizer.concrete_problem, abstract_system, optimizer.abstract_transition_cost)

    # Initial set for early stop
    init_set =
        optimizer.early_stop ? optimizer.abstract_problem.initial_set :
        SY.enum_states(abstract_system)

    optimizer.print_level >= 1 && println("compute_controller_reachability! started")

    abstract_controller,
    controllable_ids,
    uncontrollable_ids,
    value_fun_tab = SY.compute_worst_case_cost_controller(
        SY.get_automaton(abstract_system),
        optimizer.abstract_problem.target_set;
        initial_set = init_set,
        sparse_input = optimizer.sparse_input,
        cost_function = optimizer.abstract_problem.transition_cost,
    )

    optimizer.abstract_controller = abstract_controller
    optimizer.controllable_set = SY.get_state_set_from_states(abstract_system, controllable_ids)
    optimizer.uncontrollable_set = SY.get_state_set_from_states(abstract_system, uncontrollable_ids)
    optimizer.value_fun_tab = value_fun_tab

    optimizer.abstract_value_function = build_abstract_value_function(value_fun_tab)
    optimizer.concrete_value_function =
        build_concrete_value_function(abstract_system, optimizer.abstract_value_function)

    # success check: "initial_set ⊆ controllable"
    xm  = SY.get_state_mapping(abstract_system)
    optimizer.success = all(q -> MP.contains_state(optimizer.controllable_set, xm, q), init_set)

    optimizer.print_level >= 1 &&
        println("\n Reachability: terminated with $(optimizer.success)")

    optimizer.abstract_problem_time_sec = time() - t_ref
    return
end


function build_abstract_problem(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
    abstract_transition_cost
)
    @warn("The `state_cost` is not yet fully implemented")

    return Dionysos.Problem.OptimalControlProblem(
        abstract_system,
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.initial_set,
            MP.OUTER,
        ),
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.target_set,
            MP.INNER,
        ),
        concrete_problem.state_cost, # TODO: Transform continuous cost into discrete abstraction
        abstract_transition_cost,
        concrete_problem.time, # TODO: Translate continuous time into discrete steps
    )
end
