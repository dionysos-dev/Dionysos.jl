# Utils

Foundational helpers built on [`LazySets`](https://juliareach.github.io/LazySets.jl): Dionysos sets
*are* LazySets (a box is a `Hyperrectangle` via [`UT.box`](@ref Dionysos.Utils.box), an ellipsoid a
`LazySets.Ellipsoid`), plus what LazySets lacks for symbolic control — callable cost functions, exact
ellipsoid predicates ([`is_included`](@ref Dionysos.Utils.is_included) /
[`is_disjoint`](@ref Dionysos.Utils.is_disjoint)), periodic splitting, data structures,
[`RRT`](@ref Dionysos.Utils.RRT) search, and scalar optimization. The `PathCompleteFramework`
submodule provides path-complete Lyapunov functions.

## API reference

```@autodocs
Modules = [Dionysos.Utils]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Path-Complete Framework

```@autodocs
Modules = [Dionysos.Utils.PathCompleteFramework]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
