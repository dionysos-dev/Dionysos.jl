mutable struct OptimizerReachAndStayProblem{T} <: MOI.AbstractOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_system::Any
    early_stop::Bool
    print_level::Int

    # outputs
    abstract_optimizer::Union{Nothing, OPDS.OptimizerReachAndStayProblem}
    abstract_problem::Union{Nothing, PR.ReachAndStayProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    winning_set::Union{Nothing, MP.AbstractStateSet}
    winning_set_complement::Union{Nothing, MP.AbstractStateSet}
    success::Bool
    abstract_problem_time_sec::T

    function OptimizerReachAndStayProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            false,
            1,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            zero(T),
        )
    end
end

OptimizerReachAndStayProblem() = OptimizerReachAndStayProblem{Float64}()

MOI.is_empty(optimizer::OptimizerReachAndStayProblem) =
    optimizer.concrete_problem === nothing

function MOI.set(
    model::OptimizerReachAndStayProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerReachAndStayProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerReachAndStayProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstract_problem(
    concrete_problem::PR.ReachAndStayProblem,
    abstract_system::SY.SymbolicModel,
)
    return PR.ReachAndStayProblem(
        SY.get_automaton(abstract_system),
        SY.get_states_from_set(abstract_system, concrete_problem.initial_set, MP.OUTER),
        SY.get_states_from_set(abstract_system, concrete_problem.target_set, MP.INNER),
        SY.get_states_from_set(abstract_system, concrete_problem.safe_set, MP.INNER),
        concrete_problem.time,
    )
end

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing && error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
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

    optimizer.winning_set = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.winning_set),
    )

    optimizer.winning_set_complement = SY.get_state_set_from_states(
        abstract_system,
        collect(abstract_optimizer.winning_set_complement),
    )

    optimizer.success = abstract_optimizer.success
    optimizer.abstract_problem_time_sec = time() - t0

    return
end
