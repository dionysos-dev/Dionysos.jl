# Mapping

Handles **discretization** — the link between a concrete continuous state space and its finite
abstract counterpart. A [`GridFree`](@ref Dionysos.Mapping.GridFree) grid (centre `x0`, step `hx`)
defines cells, an [`AbstractMapping`](@ref Dionysos.Mapping.AbstractMapping) is the bijection between
integer state labels and cells, and cells are collected into state sets with an inclusion mode
`INNER` / `OUTER` / `CENTER`. [`MappedStateSet`](@ref Dionysos.Mapping.MappedStateSet) is the public
state-set surface; the concrete mappings trade memory for lookup speed.

## API reference

```@autodocs
Modules = [Dionysos.Mapping]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
