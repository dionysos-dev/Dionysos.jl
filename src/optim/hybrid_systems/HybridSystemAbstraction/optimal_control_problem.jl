"""
    OptimizerOptimalControlProblem{T} <: Dionysos.Optim.AbstractLiftedControlOptimizer

Reach(-avoid) sub-solver of the hybrid family: lifts an
[`OptimalControlProblem`](@ref Dionysos.Problem.OptimalControlProblem) over a hybrid system onto
the [`HybridSymbolicModel`](@ref Dionysos.Symbolic.HybridSymbolicModel) abstraction and solves it
with the discrete
[`OptimizerOptimalControlProblem`](@ref Dionysos.Optim.DiscreteSystems.OptimizerOptimalControlProblem).

A transition cost given over concrete augmented states is translated by
[`get_abstract_transition_cost`](@ref); mode switches reach it as their switching label.

Set `"concrete_problem"` and `"abstract_system"`; read back `"abstract_controller"`,
`"controllable_set"`, `"uncontrollable_set"` and `"abstract_value_function"`.
"""
mutable struct OptimizerOptimalControlProblem{T} <: AbstractLiftedControlOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, SY.HybridSymbolicModel}
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

abstract_optimizer_type(::OptimizerOptimalControlProblem) =
    OPDS.OptimizerOptimalControlProblem

function configure_abstract_optimizer!(
    model::OptimizerOptimalControlProblem,
    abstract_optimizer,
)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("early_stop"), model.early_stop)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        model.sparse_input,
    )
    return
end

function extract_results!(model::OptimizerOptimalControlProblem, abstract_optimizer)
    model.controllable_set = abstract_optimizer.controllable_set
    model.uncontrollable_set = abstract_optimizer.uncontrollable_set
    model.value_fun_tab = abstract_optimizer.value_fun_tab
    model.abstract_value_function = abstract_optimizer.value_function

    if model.print_level >= 1 && !isempty(model.abstract_problem.initial_set)
        q0 = model.abstract_problem.initial_set[1]
        println("value init_state : ", model.value_fun_tab[q0])
    end
    return
end

"""
    get_abstract_transition_cost(abstract_system, concrete_transition_cost)

Lift a transition cost written over concrete augmented states onto the abstraction, or `nothing`
if the problem carries no cost.

The returned closure is called by the discrete solver with `(abstract_state, abstract_input)` and
concretizes both before delegating: a continuous input becomes the input of the state's own mode,
while a mode switch becomes its switching label, so `cost((x, mode), "SWITCH 1 -> 2")` is how a
model prices changing mode.
"""
function get_abstract_transition_cost(
    abstract_system::SY.HybridSymbolicModel,
    concrete_transition_cost,
)
    concrete_transition_cost === nothing && return nothing

    function abstract_transition_cost(state, input)
        aug_concrete_state = SY.get_concrete_state(abstract_system, state)
        k = aug_concrete_state[end]
        if SY.is_switching_input(abstract_system.input_mapping, input)
            transition_id = abstract_system.input_mapping.global_to_switching[input]
            u = SY.get_switch_label(SY.get_global_input_map(abstract_system), transition_id)
        else
            u = SY.get_concrete_input(abstract_system, input, k)
        end
        return concrete_transition_cost(aug_concrete_state, u)
    end

    return abstract_transition_cost
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.HybridSymbolicModel,
)
    abstract_initial_set =
        _abstract_initial_states(abstract_system, concrete_problem.initial_set)
    _check_initial_nonempty(abstract_initial_set)

    abstract_target_set = SY.states_satisfying(abstract_system, concrete_problem.target_set)
    _check_nonempty(abstract_target_set, "target")

    concrete_safe_set = concrete_problem.safe_set
    abstract_safe_set =
        concrete_safe_set === nothing ? nothing :
        SY.states_satisfying(abstract_system, concrete_safe_set)
    abstract_safe_set === nothing || _check_nonempty(abstract_safe_set, "safe")

    return PR.OptimalControlProblem(
        SY.get_automaton(abstract_system),
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost, # TODO: Transform continuous cost into discrete abstraction
        get_abstract_transition_cost(abstract_system, concrete_problem.transition_cost),
        concrete_problem.time; # TODO: Translate continuous time into discrete steps
        safe_set = abstract_safe_set,
    )
end

"""
    reached(problem, aug_state) -> Bool

Whether the augmented state `(x[, t], mode)` has arrived in `problem`'s target. Both reach and
reach-and-stay ask the same question of the same field, so they share the method.
"""
reached(problem::Union{PR.OptimalControlProblem, PR.ReachAndStayProblem}, aug_state) =
    PR.satisfies(problem.target_set, aug_state...)
