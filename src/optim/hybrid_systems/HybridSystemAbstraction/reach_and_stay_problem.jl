"""
    OptimizerReachAndStayProblem{T} <: Dionysos.Optim.AbstractLiftedControlOptimizer

Reach-and-stay sub-solver of the hybrid family: lifts a
[`ReachAndStayProblem`](@ref Dionysos.Problem.ReachAndStayProblem) over a hybrid system onto the
[`HybridSymbolicModel`](@ref Dionysos.Symbolic.HybridSymbolicModel) abstraction and solves it with
the discrete
[`OptimizerReachAndStayProblem`](@ref Dionysos.Optim.DiscreteSystems.OptimizerReachAndStayProblem).

Both the target and the safe set are mode-indexed
([`HybridSpec`](@ref Dionysos.Problem.HybridSpec)), so "settle in the target" may mean a different
region in each mode. `stay_on_first_entry` carries through from the concrete problem and picks
which reading of *and stay* is enforced.

Set `"concrete_problem"` and `"abstract_system"`; read back `"abstract_controller"`,
`"winning_set"` and `"winning_set_complement"`.
"""
mutable struct OptimizerReachAndStayProblem{T} <: AbstractLiftedControlOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_system::Union{Nothing, SY.HybridSymbolicModel}
    early_stop::Bool
    print_level::Int

    # Outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerReachAndStayProblem}
    abstract_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    winning_set::Any
    winning_set_complement::Any
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerReachAndStayProblem{T}() where {T}
        return new{T}(
            nothing, # concrete_problem
            nothing, # abstract_system
            false,   # early_stop
            1,       # print_level
            nothing, # abstract_optimizer
            nothing, # abstract_problem
            nothing, # abstract_controller
            nothing, # winning_set
            nothing, # winning_set_complement
            false,   # success
            zero(T), # abstract_problem_time_sec
        )
    end
end

OptimizerReachAndStayProblem() = OptimizerReachAndStayProblem{Float64}()

function reset!(model::OptimizerReachAndStayProblem)
    model.abstract_optimizer = nothing
    model.abstract_problem = nothing
    model.abstract_controller = nothing
    model.winning_set = nothing
    model.winning_set_complement = nothing
    model.abstract_problem_time_sec = 0.0
    model.success = false
    return model
end

abstract_optimizer_type(::OptimizerReachAndStayProblem) = OPDS.OptimizerReachAndStayProblem

function configure_abstract_optimizer!(
    model::OptimizerReachAndStayProblem,
    abstract_optimizer,
)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("early_stop"), model.early_stop)
    return
end

function extract_results!(model::OptimizerReachAndStayProblem, abstract_optimizer)
    model.winning_set = abstract_optimizer.winning_set
    model.winning_set_complement = abstract_optimizer.winning_set_complement
    model.print_level >= 1 &&
        model.winning_set !== nothing &&
        println("Winning set size: $(length(model.winning_set))")
    return
end

function build_abstract_problem(
    concrete_problem::PR.ReachAndStayProblem,
    abstract_system::SY.HybridSymbolicModel,
)
    abstract_initial_set =
        _abstract_initial_states(abstract_system, concrete_problem.initial_set)
    abstract_target_set = SY.states_satisfying(abstract_system, concrete_problem.target_set)
    abstract_safe_set = SY.states_satisfying(abstract_system, concrete_problem.safe_set)

    _check_initial_nonempty(abstract_initial_set)
    _check_nonempty(abstract_target_set, "target")
    _check_nonempty(abstract_safe_set, "safe")

    return PR.ReachAndStayProblem(
        SY.get_automaton(abstract_system),
        abstract_initial_set,
        abstract_target_set,
        abstract_safe_set,
        concrete_problem.time;
        stay_on_first_entry = concrete_problem.stay_on_first_entry,
    )
end
