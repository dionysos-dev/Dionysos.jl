# ----------------------------
# Problem interface
# ----------------------------

"""
    ProblemType

Root of the control-task specification hierarchy. A problem bundles a system
with a specification; a solver consumes it through the MOI interface.

There are two categories:

- [`ControlProblem`](@ref): a control objective is synthesized (there is an
  initial set and a [`trajectory_success`](@ref) predicate).
- [`AbstractionProblem`](@ref): no control objective; the problem only
  parametrizes the construction of a (reusable) abstraction.

# Extending

To add a problem type, subtype `ControlProblem` or `AbstractionProblem` and
implement:

- [`discretize_problem`](@ref) — time-discretize the (continuous) system and
  horizon (unless the generic method already covers it via `remake`);
- [`trajectory_success`](@ref) — *control problems only* — whether a closed-loop
  trajectory satisfies the specification;
- a plotting `@recipe` (optional).
"""
abstract type ProblemType end

"""
    ControlProblem <: ProblemType

A specification with a control objective to synthesize: reach-avoid, safety,
reach-and-stay, co-safe LTL. Every `ControlProblem` has an `initial_set` and a
[`trajectory_success`](@ref) predicate.
"""
abstract type ControlProblem <: ProblemType end

"""
    AbstractionProblem <: ProblemType

A specification that carries no control objective; it parametrizes the
construction of a sound abstraction (symbolic model) that other solvers reuse.
[`trajectory_success`](@ref) is not defined for these problems.
"""
abstract type AbstractionProblem <: ProblemType end

"""
    remake(problem::ProblemType; kwargs...)

Return a copy of `problem` with the given fields replaced. Used by
[`discretize_problem`](@ref) to swap the system and horizon without a per-type
constructor call. Fields not named in `kwargs` are copied verbatim; the type
parameters are re-inferred by the constructor, so replacing a field with a value
of a different type (e.g. a discretized system) is fine.
"""
function remake(problem::P; kwargs...) where {P <: ProblemType}
    updated = NamedTuple(kwargs)
    args = ntuple(fieldcount(P)) do i
        f = fieldname(P, i)
        return haskey(updated, f) ? updated[f] : getfield(problem, f)
    end
    constructor = getfield(parentmodule(P), nameof(P))
    return constructor(args...)
end

"""
    horizon_round_up(problem::ControlProblem) -> Bool

Rounding direction used when discretizing the problem horizon: `true` for
"for at least `T`" specifications (safety, reach-and-stay), `false` for
"within at most `T`" specifications (reach-avoid). Defaults to `true`.
"""
horizon_round_up(::ControlProblem) = true

"""
    discretize_problem(problem::ProblemType, Δt::Real; num_substeps = ST.DEFAULT_NUM_SUBSTEPS)

Time-discretize `problem`: replace its continuous-time system by the
`Δt`-sampled discrete-time system and convert the horizon to a number of steps.
"""
function discretize_problem(
    problem::ProblemType,
    Δt::Real;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    return error("implement `discretize_problem` for $(typeof(problem))")
end

# Generic discretization for the control problems that only need the system and
# horizon swapped; the horizon rounding direction is the sole per-type variation.
function discretize_problem(
    problem::ControlProblem,
    Δt::Real;
    num_substeps::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    discrete_system =
        ST.discretize_continuous_system(problem.system, Δt; num_substeps = num_substeps)
    discrete_time = discretize_time(problem.time, Δt; round_up = horizon_round_up(problem))
    return remake(problem; system = discrete_system, time = discrete_time)
end

"""
    trajectory_success(problem::ControlProblem, traj::ST.Trajectory) -> Bool

Whether the closed-loop trajectory `traj` satisfies the specification of
`problem`. Defined only for [`ControlProblem`](@ref)s.
"""
trajectory_success(problem::ControlProblem, ::ST.Trajectory) =
    error("implement `trajectory_success` for $(typeof(problem))")
