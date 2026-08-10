# Utils

Foundational helpers built on [`LazySets`](https://juliareach.github.io/LazySets.jl): Dionysos sets
*are* LazySets (a box is a `Hyperrectangle`, an ellipsoid a `LazySets.Ellipsoid`), plus what LazySets
lacks for symbolic control — callable cost functions, exact ellipsoid predicates
([`is_included`](@ref Dionysos.Utils.is_included) / [`is_disjoint`](@ref Dionysos.Utils.is_disjoint)),
periodic splitting, data structures, [`RRT`](@ref Dionysos.Utils.RRT) search, and scalar optimization.
The `PathCompleteFramework` submodule provides path-complete Lyapunov functions.

## `UT.box` is a shorthand, not a set type

The examples write [`UT.box(lb, ub)`](@ref Dionysos.Utils.box) everywhere, which can read as though
Dionysos had a box type of its own. It does not. `UT.box` returns a plain
`LazySets.Hyperrectangle`; it is `LazySets.Hyperrectangle(; low, high)` with the bounds normalised to
`SVector`, and it exists mainly because the *positional* `Hyperrectangle(c, r)` takes a **centre and
a radius**, so passing bounds to it silently builds a different box.

```@example utils_box
import Dionysos
const UT = Dionysos.Utils
const LS = Dionysos.LazySets

a = UT.box([-1.0, -2.0], [1.0, 3.0])
b = LS.Hyperrectangle(; low = [-1.0, -2.0], high = [1.0, 3.0])

(LS.low(a), LS.high(a)) == (LS.low(b), LS.high(b))
```

So anywhere Dionysos asks for a set — a target, a safe set, an obstacle, a guard — you may hand it
any bounded `LazySet` and skip `UT.box` entirely:

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
