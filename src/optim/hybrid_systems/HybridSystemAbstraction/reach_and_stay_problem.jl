"""
    OptimizerReachAndStayProblem{T} <: Dionysos.Optim.AbstractDionysosOptimizer

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
mutable struct OptimizerReachAndStayProblem{T} <: OP.AbstractDionysosOptimizer
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

MOI.is_empty(optimizer::OptimizerReachAndStayProblem) =
    optimizer.concrete_problem === nothing

function MOI.get(model::OptimizerReachAndStayProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

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

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    optimizer.abstract_system === nothing && error("abstract_system not set")
    optimizer.concrete_problem === nothing && error("concrete_problem not set")

    abstract_problem =
        build_abstract_problem(optimizer.concrete_problem, optimizer.abstract_system)
    optimizer.abstract_problem = abstract_problem

    abstract_optimizer = MOI.instantiate(OPDS.OptimizerReachAndStayProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("early_stop"),
        optimizer.early_stop,
    )
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(abstract_optimizer)

    optimizer.abstract_optimizer = abstract_optimizer
    optimizer.abstract_controller = abstract_optimizer.controller
    optimizer.winning_set = abstract_optimizer.winning_set
    optimizer.winning_set_complement = abstract_optimizer.winning_set_complement
    optimizer.success = abstract_optimizer.success

    optimizer.print_level >= 1 &&
        optimizer.winning_set !== nothing &&
        println("Winning set size: $(length(optimizer.winning_set))")

    optimizer.abstract_problem_time_sec = time() - t0
    return
end

function build_abstract_problem(
    concrete_problem::PR.ReachAndStayProblem,
    abstract_system::SY.HybridSymbolicModel,
)
    q0 = SY.get_abstract_state(abstract_system, concrete_problem.initial_set)
    q0 <= 0 && error("Initial augmented state is outside the abstract system domain.")

    abstract_target_set = SY.states_satisfying(abstract_system, concrete_problem.target_set)
    abstract_safe_set = SY.states_satisfying(abstract_system, concrete_problem.safe_set)

    _check_nonempty(abstract_target_set, "target")
    _check_nonempty(abstract_safe_set, "safe")

    return PR.ReachAndStayProblem(
        SY.get_automaton(abstract_system),
        [q0],
        abstract_target_set,
        abstract_safe_set,
        concrete_problem.time;
        stay_on_first_entry = concrete_problem.stay_on_first_entry,
    )
end

"Whether the augmented state has arrived in the target — the stopping rule of a `◇□` run."
function reached(concrete_problem::PR.ReachAndStayProblem, aug_state)
    return PR.satisfies(concrete_problem.target_set, aug_state...)
end
