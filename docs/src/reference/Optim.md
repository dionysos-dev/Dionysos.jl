# Optim

The **solver catalog**. Every solver is a
[`MathOptInterface`](https://jump.dev/MathOptInterface.jl) optimizer (a subtype of
[`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer)
implementing [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!)),
so a control task is re-solved by swapping the optimizer. Solvers compose: a high-level optimizer
holds sub-solvers and forwards attributes to them.

The canonical entry point is a JuMP model, `Model(Dionysos.Optimizer)`; its vocabulary — the
operators, specification markers and mode macros — is documented in the [`Wrapper`](@ref Wrapper)
reference. Module-level entry points that belong to no submodule are in the
[`Dionysos`](Dionysos.md) reference.

## Shared solver base

Attribute forwarding, sub-solver management, and the abstraction-based composite template (compute the
abstraction, run a control sub-solver, concretize the controller).

```@autodocs
Modules = [Dionysos.Optim]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Continuous-system abstraction solvers

```@autodocs
Modules = [Dionysos.Optim.Abstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Uniform grid abstraction

SCOTS-style abstraction on a uniform grid (`GROWTH` / `LINEARIZED`), one control optimizer per
specification.

```@autodocs
Modules = [Dionysos.Optim.Abstraction.UniformGridAbstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Uniform ellipsoid abstraction

```@autodocs
Modules = [Dionysos.Optim.Abstraction.UniformEllipsoidAbstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Lazy ellipsoids abstraction

Lazy, controller-driven ellipsoidal abstraction (RRT exploration + SDP/Lyapunov transitions).

```@autodocs
Modules = [Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Hybrid-system solvers

### Hybrid system abstraction

```@autodocs
Modules = [Dionysos.Optim.Abstraction.HybridSystemAbstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### PCLF bisimulation quotient

Bisimulation-quotient synthesis for switched systems via a path-complete Lyapunov function, plus
co-safe LTL control on the quotient.

```@autodocs
Modules = [Dionysos.Optim.Abstraction.PCLFBisimulationQuotient]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Trajectory generators and certifiers

Generators produce a candidate trajectory; certifiers build a formally-certified tube around it.

```@autodocs
Modules = [Dionysos.Optim.Abstraction.OptimizerTrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

```@autodocs
Modules = [Dionysos.Optim.Abstraction.MPPITrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

```@autodocs
Modules = [Dionysos.Optim.Abstraction.CompositeTrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

```@autodocs
Modules = [Dionysos.Optim.Abstraction.UniformGridTrajectoryCertifier]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

```@autodocs
Modules = [Dionysos.Optim.Abstraction.EllipsoidalTrajectoryCertifier]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

```@autodocs
Modules = [Dionysos.Optim.Abstraction.TrajectoryCertificationOptimizer]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Discrete-system solvers

Controller synthesis directly on a finite automaton (no abstraction build).

```@autodocs
Modules = [Dionysos.Optim.DiscreteSystems]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Other solvers

### Bemporad–Morari (MIQP for PWA systems)

```@autodocs
Modules = [Dionysos.Optim.BemporadMorari]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Branch and bound

```@autodocs
Modules = [Dionysos.Optim.BranchAndBound]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
