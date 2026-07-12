# Utils

The `Utils` module is the foundational layer, built on top of
[`LazySets`](https://juliareach.github.io/LazySets.jl). Sets in Dionysos *are* LazySets: a box is a
`LazySets.Hyperrectangle` (built by [`UT.box`](@ref Dionysos.Utils.box)), an ellipsoid is a
`LazySets.Ellipsoid`. On top of that, `Utils` supplies what LazySets lacks for symbolic control:

- **Callable cost functions** — [`ScalarFunction`](@ref Dionysos.Utils.ScalarFunction),
  [`QuadraticFunction`](@ref Dionysos.Utils.QuadraticFunction),
  [`PolyhedralFunction`](@ref Dionysos.Utils.PolyhedralFunction), and black-box wrappers.
- **Geometric helpers** — the quadratic-form bridge
  ([`get_quadratic_form`](@ref Dionysos.Utils.get_quadratic_form)), sublevel sets, sampling, and
  periodic splitting ([`set_in_period`](@ref Dionysos.Utils.set_in_period)).
- **Set algebra** — unions, set-minus, and the Dionysos-owned set predicates
  [`is_included`](@ref Dionysos.Utils.is_included) / [`is_disjoint`](@ref Dionysos.Utils.is_disjoint)
  (used instead of piracy on `Base` methods over LazySets types).
- **Data structures & search** — trees, sorted vector sets, and
  [`RRT`](@ref Dionysos.Utils.RRT).
- **Scalar optimization** — [`golden_section_search`](@ref Dionysos.Utils.golden_section_search),
  [`newton_method`](@ref Dionysos.Utils.newton_method),
  [`derivative_bisection`](@ref Dionysos.Utils.derivative_bisection).
- **Path-Complete Framework** ([`PathCompleteFramework`](@ref Dionysos.Utils.PathCompleteFramework)) —
  path-complete Lyapunov functions over labelled digraphs, used by the PCLF bisimulation-quotient
  solver.

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
