# Problem

The `Problem` module defines the **specifications** — the control tasks Dionysos can solve. A
specification is deliberately *solver-independent*: the same problem can be handed to any compatible
solver in [`Optim`](@ref Optim), which is what makes algorithms swappable and comparable.

Every specification is a subtype of [`ProblemType`](@ref Dionysos.Problem.ProblemType), which splits
into two families:

- [`ControlProblem`](@ref Dionysos.Problem.ControlProblem) — a controller is synthesized to meet a
  behavioural objective. These carry an initial set and a
  [`trajectory_success`](@ref Dionysos.Problem.trajectory_success) predicate that decides whether a
  closed-loop run satisfies the spec:
  - [`OptimalControlProblem`](@ref Dionysos.Problem.OptimalControlProblem) — reach-avoid over a finite
    horizon, with optional state and transition costs.
  - [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem) — stay inside a safe set for the whole
    horizon.
  - [`ReachAndStayProblem`](@ref Dionysos.Problem.ReachAndStayProblem) — eventually reach and then
    remain in a target set.
  - [`CoSafeLTLProblem`](@ref Dionysos.Problem.CoSafeLTLProblem) — satisfy a co-safe LTL formula
    (reach an accepting condition in finite time, i.e. achieve a "good prefix").
- [`AbstractionProblem`](@ref Dionysos.Problem.AbstractionProblem) — no control objective; the problem
  only parametrizes the construction of a reusable abstraction:
  - [`AlternatingSimulationProblem`](@ref Dionysos.Problem.AlternatingSimulationProblem) — build a
    sound abstraction of a system.
  - [`BisimulationQuotientProblem`](@ref Dionysos.Problem.BisimulationQuotientProblem) — build a
    quotient bisimulation of a switched system from an observation map.

Infinite horizons are expressed with the [`Infinity`](@ref Dionysos.Problem.Infinity) sentinel.

## API reference

```@autodocs
Modules = [Dionysos.Problem]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
