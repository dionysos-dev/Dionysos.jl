# Utils

Foundational helpers built on [`LazySets`](https://juliareach.github.io/LazySets.jl): Dionysos sets
*are* LazySets (a box is a `Hyperrectangle`, an ellipsoid a `LazySets.Ellipsoid`), plus what LazySets
lacks for symbolic control — callable cost functions, exact ellipsoid predicates
([`is_included`](@ref Dionysos.Utils.is_included) / [`is_disjoint`](@ref Dionysos.Utils.is_disjoint)),
periodic splitting, data structures, [`RRT`](@ref Dionysos.Utils.RRT) search, and scalar optimization.
The `PathCompleteFramework` submodule provides path-complete Lyapunov functions.

## Sets are LazySets, with no Dionysos wrapper

There is no Dionysos box type and no Dionysos box constructor. A box is a
`LazySets.Hyperrectangle`, written from its bounds:

```@example utils_box
import Dionysos
const LS = Dionysos.LazySets

LS.Hyperrectangle(; low = [-1.0, -2.0], high = [1.0, 3.0])
```

Always the **keyword** form. Positionally, `Hyperrectangle(c, r)` takes a *centre and a radius*, so
the same numbers build a different box and nothing warns you:

```@example utils_box
box = LS.Hyperrectangle(; low = [-1.0, -2.0], high = [1.0, 3.0])
trap = LS.Hyperrectangle([-1.0, -2.0], [1.0, 3.0])

(LS.low(box), LS.high(box)), (LS.low(trap), LS.high(trap))
```

Anywhere Dionysos asks for a set — a target, a safe set, an obstacle, a guard — any bounded
`LazySet` will do:

```@example utils_box
LS.BallInf([0.0, 0.0], 1.0), LS.Ball2([0.0, 0.0], 1.0), LS.Zonotope([1.0, 1.0], [0.3 0.1; 0.0 0.2])
```

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
