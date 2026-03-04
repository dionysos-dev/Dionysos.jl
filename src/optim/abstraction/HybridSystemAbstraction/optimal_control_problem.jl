mutable struct OptimizerOptimalControlProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, SY.TimedHybridSymbolicModel}

    # Specific parameters
    optimizer_kwargs_dict::Union{Nothing, Any}

    print_level::Int

    # Outputs
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_controller::Union{Nothing, MS.AbstractMap}
    abstract_problem_time_sec::T
    success::Bool

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            1,
            nothing,
            nothing,
            0.0,
            false
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
    model.abstract_controller = nothing
    model.abstract_problem_time_sec = 0.0
    model.success = false
    return model
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t_ref = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing &&
        error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    # Build abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem


    optimizer.print_level >= 1 && println("compute_controller_reachability! started")

    abstract_controller, controllable_set_symbols, _, value_per_node = SY.compute_worst_case_cost_controller(
                abstract_problem.system.symbolic_automaton,
                abstract_problem.target_set;
                initial_set = abstract_problem.initial_set,
                sparse_input = false,
                cost_function = abstract_problem.transition_cost,
    )

    optimizer.abstract_controller = abstract_controller
    optimizer.success = ⊆(abstract_problem.initial_set, controllable_set_symbols)

    if optimizer.success
        optimizer.print_level >= 1 && println("✅ Optimal control problem is solvable: initial set is controllable")
    else
        optimizer.print_level >= 1 && println("⚠️ Warning: initial set is only partially controllable")
    end
    optimizer.print_level >= 1 && println("value init_state : ", value_per_node[abstract_problem.initial_set[1]])
    optimizer.abstract_problem_time_sec = time() - t_ref
    return
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.TimedHybridSymbolicModel,
)
    concrete_initial_state = concrete_problem.initial_set # a unique augmented point
    abstract_initial_set = [
        SY.get_abstract_state(
            abstract_system,
            concrete_initial_state,
        ),
    ]

    concrete_target_set = concrete_problem.target_set
    abstract_target_set = SY.get_states_from_set(
        abstract_system,
        concrete_target_set[1], # state
        concrete_target_set[2], # time
        concrete_target_set[3], # mode
    )

    return PR.OptimalControlProblem(
        abstract_system,
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost, # TODO: Transform continuous cost into discrete abstraction
        get_abstract_transition_cost(abstract_system, concrete_problem.transition_cost),
        concrete_problem.time, # TODO: Translate continuous time into discrete steps
    )
end

function get_abstract_transition_cost(
    abstract_system::SY.TimedHybridSymbolicModel,
    concrete_transition_cost,
)
    function abstract_transition_cost(state, input)
        (x, t, k) = SY.get_concrete_state(
            abstract_system,
            state,
        )
        aug_concrete_state = (x, t, k)
        if SY.is_switching_input(
            abstract_system.input_mapping,
            input,
        )
            transition_id = abstract_system.input_mapping.global_to_switching[input]
            label = abstract_system.input_mapping.switch_labels[transition_id]
            u = label
        else
            u = SY.get_concrete_input(
                abstract_system,
                input,
                k,
            )
        end
        return concrete_transition_cost(aug_concrete_state, u)
    end
    return abstract_transition_cost
end

function reached(concrete_problem::PR.OptimalControlProblem, aug_state)
    (x, t, k) = aug_state
    (Xs_target, Ts_target, Ns_target) = concrete_problem.target_set
    idx = findfirst(==(k), Ns_target)
    if isnothing(idx)
        return false
    end
    X_set = Xs_target[idx]
    T_set = Ts_target[idx]
    in_X = x ∈ X_set
    in_T = t ≥ T_set.lb[1] && t ≤ T_set.ub[1]

    return in_X && in_T
end
