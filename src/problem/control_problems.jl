# ----------------------------
# Control problems
# ----------------------------

"""
    OptimalControlProblem{S, XI, XT, XC, TC, T} <: ControlProblem

Encodes a **reach-avoid optimal control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached.
- `XC`: A state cost function or structure.
- `TC`: A transition cost function or structure.
- `T`: Satisfy the property in at most time T.

This problem aims to find a control strategy that reaches the target set from the initial set, minimizing the accumulated cost over time.
"""
struct OptimalControlProblem{S, XI, XT, XC, TC, T} <: ControlProblem
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
end

OptimalControlProblem(system, initial_set, target_set, state_cost, transition_cost) =
    OptimalControlProblem(
        system,
        initial_set,
        target_set,
        state_cost,
        transition_cost,
        Infinity(),
    )

# Reach-avoid is a "within at most T" specification: round the horizon down.
horizon_round_up(::OptimalControlProblem) = false

function trajectory_success(problem::OptimalControlProblem, traj::ST.Trajectory)
    isempty(traj.states) && return false

    return first(traj.states) ∈ problem.initial_set &&
           any(x -> x ∈ problem.target_set, traj.states)
end

"""
    SafetyProblem{S, XI, XS, T} <: ControlProblem

Encodes a **safety control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XS`: The safe set in which the system must remain.
- `T`: Satisfy the property for at least time T.

This problem aims to synthesize a controller that ensures the system remains within the safe set for the entire duration of the time horizon.
"""
struct SafetyProblem{S, XI, XS, T} <: ControlProblem
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

SafetyProblem(system, initial_set, safe_set) =
    SafetyProblem(system, initial_set, safe_set, Infinity())

function trajectory_success(problem::SafetyProblem, traj::ST.Trajectory)
    isempty(traj.states) && return false

    return first(traj.states) ∈ problem.initial_set &&
           all(x -> x ∈ problem.safe_set, traj.states)
end

"""
    ReachAndStayProblem{S, XI, XT, XS, T} <: ControlProblem

Encodes a **reach-and-stay control problem** (eventually always).

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached and stayed in.
- `XS`: The safe set in which the system must remain during the approach.
- `T`: Satisfy the property for at least time T.

This problem aims to synthesize a controller that drives the system from the initial set
into the target set and keeps it there indefinitely, while remaining within the safe set
during the approach phase.
"""
struct ReachAndStayProblem{S, XI, XT, XS, T} <: ControlProblem
    system::S
    initial_set::XI
    target_set::XT
    safe_set::XS
    time::T
end

ReachAndStayProblem(system, initial_set, target_set, safe_set) =
    ReachAndStayProblem(system, initial_set, target_set, safe_set, Infinity())

function trajectory_success(problem::ReachAndStayProblem, traj::ST.Trajectory)
    xs = traj.states
    isempty(xs) && return false

    first(xs) ∈ problem.initial_set || return false
    all(x -> x ∈ problem.safe_set, xs) || return false

    # Eventually-always-in-target: some suffix `xs[k:end]` lies entirely in the
    # target. Every such suffix contains `xs[end]`, and the length-1 suffix at
    # `k = length(xs)` witnesses it, so this reduces to the last state being in
    # the target.
    return last(xs) ∈ problem.target_set
end

"""
    CoSafeLTLProblem{S, XI, SPEC, LAB} <: ControlProblem

Encodes a **co-safe LTL control problem**.

- `S`: The system to control.
- `XI`: The initial set of states.
- `SPEC`: The co-safe LTL specification object (an automaton/monitor wrapper).
- `LAB`: The labeling payload type used in `labeling` (typically a concrete set type such as a LazySet,
         or an abstract labeling such as `Vector{Int}` / bitset / etc.).

# Fields
- `system::S`:
  The (concrete or abstract) system to control.

- `initial_set::XI`:
  Initial set of states (or initial abstract states).

- `spec::SPEC`:
  The co-safe LTL specification.

- `labeling::Dict{Symbol, LAB}`:
  Unified container mapping each atomic proposition (AP) `:ap` to its labeling object.
  In a **concrete** problem, values are typically sets (e.g. LazySets / Dionysos sets) over the state space.
  In an **abstract** problem, values are typically collections of abstract states (e.g. `Vector{Int}`).

- `ap_semantics::Dict{Symbol, UT.INCL_MODE}`:
  Per-AP semantics used when converting set labels to abstract labels
  (`UT.INNER` or `UT.OUTER`; also reachable as `Mapping.INNER`/`Mapping.OUTER`).


This problem aims to synthesize a controller such that the generated trajectory satisfies the co-safe LTL
formula, i.e. it reaches an accepting condition in finite time.
"""
struct CoSafeLTLProblem{S, XI, SPEC, LAB} <: ControlProblem
    system::S
    initial_set::XI
    spec::SPEC

    # unified labeling container:
    labeling::Dict{Symbol, LAB}   # Symbol => LazySet (concrete) or Vector{Int} (abstract)

    ap_semantics::Dict{Symbol, UT.INCL_MODE}
end

# Parametric structs do not convert non-parametric fields, so accept any
# INCL_MODE-valued dict (e.g. the common `Dict{Symbol, Any}`) and convert.
function CoSafeLTLProblem(
    system,
    initial_set,
    spec,
    labeling::Dict{Symbol, LAB},
    ap_semantics::AbstractDict{Symbol},
) where {LAB}
    return CoSafeLTLProblem(
        system,
        initial_set,
        spec,
        labeling,
        convert(Dict{Symbol, UT.INCL_MODE}, ap_semantics),
    )
end

function trajectory_success(problem::CoSafeLTLProblem, traj::ST.Trajectory)
    isempty(traj.states) && return false

    # Placeholder until monitor/spec trajectory evaluation is implemented.
    return false
end
