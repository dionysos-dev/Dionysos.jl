# ----------------------------------------------------------------------------------------
# Results: solution status and closed-loop simulation.
#
# The status attributes are answered by the solvers themselves
# (`Optim.AbstractionControlOptimizer`), so the front-end only forwards them: a JuMP user and a
# direct-MOI user see the same status for the same model.
# ----------------------------------------------------------------------------------------

const _FORWARDED_STATUS = Union{
    MOI.TerminationStatus,
    MOI.PrimalStatus,
    MOI.DualStatus,
    MOI.ResultCount,
    MOI.RawStatusString,
    MOI.SolveTimeSec,
}

MOI.get(model::Optimizer, attr::_FORWARDED_STATUS) = MOI.get(model.inner, attr)

# ----------------------------------------------------------------------------------------
# Closed-loop simulation
# ----------------------------------------------------------------------------------------

# Stopping criterion implied by the specification: reach problems stop on reaching the
# target, safety problems stop on leaving the safe set.
_stopping_for(p::PR.OptimalControlProblem) = x -> x ∈ p.target_set
_stopping_for(p::PR.ReachAndStayProblem) = x -> x ∈ p.target_set
_stopping_for(p::PR.SafetyProblem) = x -> x ∉ p.safe_set
_stopping_for(::Any) = _ -> false

"""
    simulate(model, x0; nsteps = 100, stopping = nothing) -> System.Trajectory

Run the synthesized controller in closed loop from `x0` for at most `nsteps` steps.

`model` is the JuMP model (or the underlying [`Optimizer`](@ref)) *after* `optimize!`. The
stopping criterion defaults to the one implied by the specification — reaching the target, or
leaving the safe set — and can be overridden with `stopping`.

Returns the channelled [`Dionysos.System.Trajectory`](@ref); read it with `System.states`,
`System.inputs`, … or plot it directly.
"""
const HSA = OP.Abstraction.HybridSystemAbstraction

# A hybrid closed loop steps the augmented state `(x[, t], mode)` and needs the abstraction
# step of each mode, so it runs on its own engine.
function _simulate_hybrid(model::Optimizer, controller, aug0, nsteps, stopping)
    problem = model.problem
    stop = if stopping !== nothing
        stopping
    elseif problem isa PR.SafetyProblem
        aug -> !HSA.safe(problem, aug)
    else
        aug -> HSA.reached(problem, aug)
    end

    tsteps = [_mode_time_step(model, k) for k in mode_ids(model.ir)]
    aug_traj, u_traj = HSA.get_closed_loop_trajectory(
        problem.system,
        controller,
        tsteps,
        aug0,
        nsteps;
        stopping = stop,
    )
    return HSA.channelled_trajectory(aug_traj, u_traj)
end

function _mode_time_step(model::Optimizer, id::Int)
    for (name, value) in Iterators.reverse(model.ir.modes[id].attributes)
        name == "time_step" && return value
    end
    step = _time_step(model)
    step === nothing && error("Mode $id has no `time_step`; simulation needs one per mode.")
    return step
end

function simulate(model::Optimizer, x0; nsteps::Int = 100, stopping = nothing)
    model.solved || error("`simulate` needs a solved model: call `optimize!` first.")

    controller = MOI.get(model.inner, MOI.RawOptimizerAttribute("concrete_controller"))
    controller === nothing &&
        error("The model has no concrete controller; it carries no control objective.")

    if is_hybrid(model.ir)
        return _simulate_hybrid(model, controller, x0, nsteps, stopping)
    end

    system = MOI.get(model.inner, MOI.RawOptimizerAttribute("discrete_time_system"))
    stop = stopping === nothing ? _stopping_for(model.problem) : stopping

    return ST.get_closed_loop_trajectory(system, controller, x0, nsteps; stopping = stop)
end

function simulate(model::JuMP.GenericModel, x0; kwargs...)
    inner = JuMP.unsafe_backend(model)
    inner isa Optimizer || error(
        "`simulate` expects a model built with `Model(Dionysos.Optimizer)`, got $(typeof(inner)).",
    )
    return simulate(inner, x0; kwargs...)
end
