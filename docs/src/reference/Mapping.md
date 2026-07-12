# Mapping

The `Mapping` module handles **discretization**: it links a concrete continuous state space to its
abstract, finite counterpart and back. The space is covered by a **grid**
([`GridFree`](@ref Dionysos.Mapping.GridFree) with centre `x0` and step `hx`); a **cell** is an
integer coordinate `pos` standing for the point `x0 + hx .* pos`; and an
[`AbstractMapping`](@ref Dionysos.Mapping.AbstractMapping) is the bijection between integer state
labels `1:n` and cells. Cells are collected into state sets with an inclusion mode
`INNER` / `OUTER` / `CENTER` that fixes how a continuous region is over- or under-approximated by
cells.

The concrete mapping implementations trade memory for lookup speed:
[`ExplicitGridMapping`](@ref Dionysos.Mapping.ExplicitGridMapping) materializes the label ↔ cell
table, [`ImplicitGridMapping`](@ref Dionysos.Mapping.ImplicitGridMapping) computes it on the fly,
[`PeriodicGridMapping`](@ref Dionysos.Mapping.PeriodicGridMapping) wraps periodic coordinates, and the
[`HierarchicalGridMapping`](@ref Dionysos.Mapping.HierarchicalGridMapping) supports multi-level
(lazy) abstractions. [`MappedStateSet`](@ref Dionysos.Mapping.MappedStateSet) is the public state-set
surface — a set of cells bundled with the mapping that gives them meaning.

## API reference

```@autodocs
Modules = [Dionysos.Mapping]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
