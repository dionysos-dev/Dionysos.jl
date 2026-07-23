# Symbolic

Builds the finite **automaton abstraction** of a concrete system on top of a
[`Mapping`](@ref Mapping). A [`SymbolicModel`](@ref Dionysos.Symbolic.SymbolicModel) (concretely a
[`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList)) wraps an automaton whose transitions
are a sound over-approximation of the dynamics. The transition relation is populated by an execution
backend (sequential or parallel). Optional, composable *lifts* —
[`ClockLift`](@ref Dionysos.Symbolic.ClockLift) (time) and
[`ModeLift`](@ref Dionysos.Symbolic.ModeLift) (hybrid modes) — augment a base model into `(x, t)`,
`(x, k)`, or `(x, k, t)` abstractions.

## API reference

```@autodocs
Modules = [Dionysos.Symbolic]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
