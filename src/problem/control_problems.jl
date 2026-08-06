# ----------------------------
# Control problems
# ----------------------------

"""
    OptimalControlProblem{S, XI, XT, XC, TC, T, XS} <: ControlProblem

Encodes a **reach-avoid optimal control problem** over a finite horizon.

- `S`: The system to control.
- `XI`: The initial set of states.
- `XT`: The target set to be reached.
- `XC`: A state cost function or structure.
- `TC`: A transition cost function or structure.
- `T`: Satisfy the property in at most time T.
- `XS`: The safe set the trajectory must stay in until the target is reached, or `nothing`
  for the whole state space.

This problem aims to find a control strategy that reaches the target set from the initial set, minimizing the accumulated cost over time.

# The safe set

`safe_set` is the ◻ of a reach-avoid specification `◻ safe ∧ ◇ target`, and is optional:
`nothing` — the default, and what the five-and six-argument constructors give — means the
whole state space, so a problem built without it behaves exactly as before.

It is *not* the same as carving the unsafe region out of `stateset(system)`. A region removed
from the state space is never abstracted, so the synthesis cannot reason about it; a safe set
keeps it representable and lets the synthesis actively avoid it.

Semantics: `safe U target`. Every state up to **and including** the one that reaches the
target must be safe, so the target is effectively intersected with the safe set — a target
state outside `safe_set` does not count as reached.

```julia
PR.OptimalControlProblem(system, X0, target, nothing, nothing, PR.Infinity(); safe_set = corridor)
```
"""
struct OptimalControlProblem{S, XI, XT, XC, TC, T, XS} <: ControlProblem
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
    safe_set::XS
end

function OptimalControlProblem(
    system,
    initial_set,
    target_set,
    state_cost,
    transition_cost,
    time;
    safe_set = nothing,
)
    return OptimalControlProblem(
        system,
        initial_set,
        target_set,
        state_cost,
        transition_cost,
        time,
        safe_set,
    )
end

function OptimalControlProblem(
    system,
    initial_set,
    target_set,
    state_cost,
    transition_cost;
    safe_set = nothing,
)
    return OptimalControlProblem(
        system,
        initial_set,
        target_set,
        state_cost,
        transition_cost,
        Infinity();
        safe_set = safe_set,
    )
end

# Reach-avoid is a "within at most T" specification: round the horizon down.
horizon_round_up(::OptimalControlProblem) = false

"""
    check_safe_set_supported(problem::OptimalControlProblem, solver_name)

Raise if `problem` carries a [`safe_set`](@ref OptimalControlProblem) the caller cannot
honour.

Solvers that ignore the avoid part of a reach-avoid specification must call this. Dropping it
silently would hand back a controller certified against a weaker specification than the one
that was asked for, which is exactly the kind of failure a formal-methods toolbox must not
have.
"""
function check_safe_set_supported(problem::OptimalControlProblem, solver_name)
    problem.safe_set === nothing && return
    return error(
        "$solver_name cannot honour the `safe_set` of a reach-avoid problem, and refuses to " *
        "ignore it silently. Remove the safe set, fold it into the system's state set, or " *
        "use an abstraction-based solver.",
    )
end

function trajectory_success(problem::OptimalControlProblem, traj::ST.Trajectory)
    xs = traj.states
    isempty(xs) && return false
    first(xs) ∈ problem.initial_set || return false

    # `safe U target`: the run stops at the target, so what happens afterwards is outside the
    # specification — but every state up to and including it must be safe.
    for x in xs
        problem.safe_set === nothing || x ∈ problem.safe_set || return false
        x ∈ problem.target_set && return true
    end
    return false
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
