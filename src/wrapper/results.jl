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

# A periodic solver folds the state space into one period, so the controller is only defined
# there. The closed loop has to fold the same way — otherwise the state walks out of the
# controller's domain after one revolution — and the specification sets have to be compared in
# the same coordinates. `nothing` when the solver was not configured periodically.
struct PeriodicFold{P, T}
    dims::SVector{P, Int}
    periods::SVector{P, T}
    start::SVector{P, T}
end

(f::PeriodicFold)(x) = UT.wrap_coord(x, f.dims, f.periods; start = f.start)
_fold_set(f::PeriodicFold, set) = UT.set_in_period(set, f.dims, f.periods, f.start)
_fold_set(::Nothing, set) = set

# Attributes are field-backed: an optimizer without periodic support throws rather than
# answering `nothing`, so every read is guarded.
function _raw(model::Optimizer, name, default)
    return try
        v = MOI.get(model.inner, MOI.RawOptimizerAttribute(name))
        v === nothing ? default : v
    catch
        default
    end
end

function _periodic_fold(model::Optimizer)
    _raw(model, "use_periodic_mapping", false) === true || return nothing
    dims = _raw(model, "periodic_dims", nothing)
    periods = _raw(model, "periodic_periods", nothing)
    (dims === nothing || periods === nothing) && return nothing
    start = _raw(model, "periodic_start", zeros(SVector{length(dims), Float64}))
    return PeriodicFold(SVector(dims...), SVector(periods...), SVector(start...))
end

# Stopping criterion implied by the specification: reach problems stop on reaching the
# target, safety problems stop on leaving the safe set.
function _stopping_for(p::PR.OptimalControlProblem, fold = nothing)
    target = _fold_set(fold, p.target_set)
    p.safe_set === nothing && return x -> x ∈ target
    # Reach-avoid: the run is settled either way once the state leaves the safe set — no
    # continuation can satisfy `safe U target` — so stop there rather than simulate on.
    safe = _fold_set(fold, p.safe_set)
    return x -> x ∈ target || x ∉ safe
end
function _stopping_for(p::PR.ReachAndStayProblem, fold = nothing)
    target = _fold_set(fold, p.target_set)
    return x -> x ∈ target
end
function _stopping_for(p::PR.SafetyProblem, fold = nothing)
    safe = _fold_set(fold, p.safe_set)
    return x -> x ∉ safe
end
_stopping_for(::Any, fold = nothing) = _ -> false

# Whether the closed loop must feed the *next* measurement to the controller's memory update.
# A co-safe controller's memory is the specification automaton, and that automaton reads the
# labels of the state the system lands in. Handing it the current state instead desynchronises
# the monitor by one step, which does not fail loudly — it silently truncates the run and
# reports a trajectory that never satisfies the formula.
_update_on_next(::PR.CoSafeLTLProblem) = true
_update_on_next(::Any) = false

const HSA = OP.Abstraction.HybridSystemAbstraction

# A hybrid closed loop steps the augmented state `(x[, t], mode)` and needs the abstraction
# step of each mode, so it runs on its own engine.
function _simulate_hybrid(model::Optimizer, controller, aug0, nsteps, stopping)
    problem = model.problem
    stop = if stopping !== nothing
        stopping
    elseif problem isa PR.SafetyProblem
        aug -> !HSA.safe(problem, aug)
    elseif problem isa PR.ReachAndStayProblem
        # Same rule as the continuous path (`_stopping_for`): stop on arrival. For `◇□` that
        # hides the staying, so pass `stopping = _ -> false` to watch a run settle.
        aug -> HSA.reached(problem, aug) || !PR.satisfies(problem.safe_set, aug...)
    elseif problem isa PR.OptimalControlProblem && problem.safe_set !== nothing
        aug -> HSA.reached(problem, aug) || !PR.satisfies(problem.safe_set, aug...)
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

"""
    simulate(model, x0; nsteps = 100, stopping = nothing) -> System.Trajectory

Run the synthesized controller in closed loop from `x0` for at most `nsteps` steps.

`model` is the JuMP model (or the underlying [`Optimizer`](@ref)) *after* `optimize!`. The
stopping criterion defaults to the one implied by the specification — reaching the target, or
leaving the safe set — and can be overridden with `stopping`.

When the solver was configured with `use_periodic_mapping`, the closed loop folds the state into
the same period the abstraction used, and the specification sets are compared in those
coordinates. Without it the state would leave the controller's domain after one revolution.

Returns the channelled [`Dionysos.System.Trajectory`](@ref); read it with `System.states`,
`System.inputs`, … or plot it directly.
"""
function simulate(model::Optimizer, x0; nsteps::Int = 100, stopping = nothing)
    model.solved || error("`simulate` needs a solved model: call `optimize!` first.")

    controller = MOI.get(model.inner, MOI.RawOptimizerAttribute("concrete_controller"))
    controller === nothing &&
        error("The model has no concrete controller; it carries no control objective.")

    if is_hybrid(model.ir)
        return _simulate_hybrid(model, controller, x0, nsteps, stopping)
    end

    system = MOI.get(model.inner, MOI.RawOptimizerAttribute("discrete_time_system"))
    fold = _periodic_fold(model)
    stop = stopping === nothing ? _stopping_for(model.problem, fold) : stopping
    wrap = fold === nothing ? identity : fold

    return ST.get_closed_loop_trajectory(
        system,
        controller,
        x0,
        nsteps;
        stopping = stop,
        wrap = wrap,
        update_on_next = _update_on_next(model.problem),
    )
end

function simulate(model::JuMP.GenericModel, x0; kwargs...)
    inner = JuMP.unsafe_backend(model)
    inner isa Optimizer || error(
        "`simulate` expects a model built with `Model(Dionysos.Optimizer)`, got $(typeof(inner)).",
    )
    return simulate(inner, x0; kwargs...)
end
