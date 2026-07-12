# Symbolic

The `Symbolic` module builds the finite **automaton abstraction** of a concrete system on top of a
[`Mapping`](@ref Mapping). A [`SymbolicModel`](@ref Dionysos.Symbolic.SymbolicModel) (concretely a
[`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList)) wraps an automaton whose transitions
`(target, source, symbol)` form a sound over-approximation of the dynamics, so a controller
synthesized on it is valid on the original system.

The pieces:

- **Symbolic models** — the abstraction interface and its grid-based implementations.
- **Automaton** — the transition graph; the concrete list types trade memory for `pre` / `post`
  speed.
- **Execution backends** — how the transition relation is populated from the concrete system,
  sequentially or in parallel ([`ThreadedBackend`](@ref Dionysos.Symbolic.ThreadedBackend),
  [`JuliaDistributedBackend`](@ref Dionysos.Symbolic.JuliaDistributedBackend),
  [`SlurmArrayBackend`](@ref Dionysos.Symbolic.SlurmArrayBackend)).
- **Timed hybrid symbolic model** — abstraction of a timed hybrid system: per-mode spatial and time
  abstractions flattened into a single automaton, with inputs unified through a global input map. The
  optimizer-driven builders live in
  [`HybridSystemAbstraction`](@ref Dionysos.Optim.Abstraction.HybridSystemAbstraction); the data
  structures live here.

## API reference

```@autodocs
Modules = [Dionysos.Symbolic]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
