mutable struct OptimizerOptimalControlProblem{T} <: OP.AbstractDionysosOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, SY.TimedHybridSymbolicModel}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # Outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerOptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    controllable_set::Any
    uncontrollable_set::Any
    value_fun_tab::Any
    abstract_value_function::Any
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
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
    model.controllable_set = nothing
    model.uncontrollable_set = nothing
    model.value_fun_tab = nothing
    model.abstract_value_function = nothing
    model.abstract_problem_time_sec = 0.0
    model.success = false
    return model
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t0 = time()

    optimizer.abstract_system === nothing && error("abstract_system not set")
    optimizer.concrete_problem === nothing && error("concrete_problem not set")

    abs_sys = optimizer.abstract_system
    concrete_problem = optimizer.concrete_problem

    abstract_problem = build_abstract_problem(concrete_problem, abs_sys)
    optimizer.abstract_problem = abstract_problem

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
    optimizer.controllable_set = abstract_optimizer.controllable_set
    optimizer.uncontrollable_set = abstract_optimizer.uncontrollable_set
    optimizer.value_fun_tab = abstract_optimizer.value_fun_tab
    optimizer.abstract_value_function = abstract_optimizer.value_function
    optimizer.success = abstract_optimizer.success

    if optimizer.print_level >= 1 && !isempty(abstract_problem.initial_set)
        q0 = abstract_problem.initial_set[1]
        println("value init_state : ", optimizer.value_fun_tab[q0])
    end

    optimizer.abstract_problem_time_sec = time() - t0
    return
end

function get_abstract_transition_cost(
    abstract_system::SY.TimedHybridSymbolicModel,
    concrete_transition_cost,
)
    concrete_transition_cost === nothing && return nothing

    function abstract_transition_cost(state, input)
        (x, t, k) = SY.get_concrete_state(abstract_system, state)
        aug_concrete_state = (x, t, k)
        if SY.is_switching_input(abstract_system.input_mapping, input)
            transition_id = abstract_system.input_mapping.global_to_switching[input]
            label = abstract_system.input_mapping.switch_labels[transition_id]
            u = label
        else
            u = SY.get_concrete_input(abstract_system, input, k)
        end
        return concrete_transition_cost(aug_concrete_state, u)
    end

    return abstract_transition_cost
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.TimedHybridSymbolicModel,
)
    concrete_initial_state = concrete_problem.initial_set # a unique augmented point
    q0 = SY.get_abstract_state(abstract_system, concrete_initial_state)
    q0 === nothing &&
        error("Initial augmented state is outside the abstract system domain.")
    abstract_initial_set = [q0]

    abstract_target_set = SY.states_satisfying(abstract_system, concrete_problem.target_set)

    return PR.OptimalControlProblem(
        SY.get_automaton(abstract_system),
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost, # TODO: Transform continuous cost into discrete abstraction
        get_abstract_transition_cost(abstract_system, concrete_problem.transition_cost),
        concrete_problem.time, # TODO: Translate continuous time into discrete steps
    )
end

function reached(concrete_problem::PR.OptimalControlProblem, aug_state)
    (x, t, k) = aug_state
    return PR.satisfies(concrete_problem.target_set, x, t, k)
end
