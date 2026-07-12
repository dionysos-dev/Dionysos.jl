# Symbolic

Builds the finite **automaton abstraction** of a concrete system on top of a
[`Mapping`](@ref Mapping). A [`SymbolicModel`](@ref Dionysos.Symbolic.SymbolicModel) (concretely a
[`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList)) wraps an automaton whose transitions
are a sound over-approximation of the dynamics. The transition relation is populated by an execution
backend (sequential or parallel), and a dedicated timed-hybrid symbolic model abstracts timed hybrid
systems.

## API reference

```@autodocs
Modules = [Dionysos.Symbolic]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
