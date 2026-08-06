# Dionysos

The top-level module surface: the entry points that belong to no submodule. Everything else is
documented under [`Utils`](@ref Utils), [`System`](@ref System), [`Problem`](@ref Problem),
[`Mapping`](@ref Mapping), [`Symbolic`](@ref Symbolic), [`Optim`](@ref Optim) and
[`Wrapper`](@ref Wrapper).

Several of these are **stubs backed by package extensions**: the method only exists once the
optional dependency is loaded, and calling it before that raises an error saying so.

## Visualisation

```@docs
Dionysos.animate_trajectory_dashboard
Dionysos.animate_mechanism_trajectory
Dionysos.plot_lifted_bisimulation!
Dionysos.plot_lifted_trajectory!
```

## Import and export

```@docs
Dionysos.export_controller_csv
Dionysos.import_controller_csv
```

## Temporal logic

```@docs
Dionysos.spot_stepper
```
