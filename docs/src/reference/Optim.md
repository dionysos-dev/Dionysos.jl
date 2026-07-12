# Optim

The `Optim` module is the **solver catalog**. Every solver — abstraction-based or not — is a
[`MathOptInterface`](https://jump.dev/MathOptInterface.jl) optimizer (a subtype of
[`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer)
implementing [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!)).
This shared contract is the architectural keystone: a control task can be re-solved, compared, and
benchmarked by *swapping the optimizer* rather than rewriting the model. Solvers also **compose** — a
high-level optimizer holds sub-solvers (an abstraction solver and a control solver) and forwards
attribute `set` / `get` to them.

## JuMP frontend

The canonical user entry point is a JuMP model, `Model(Dionysos.Optimizer)`. Dynamics and target /
initial constraints are written with the `NonlinearOperator`s below; `Dionysos.Optimizer` dispatches
the model to a concrete solver family. (`Dionysos.Optimizer` is enabled by the
`DionysosMathOptSymbolicAD` extension.)

```@autodocs
Modules = [Dionysos]
Filter  = _is_public
Order   = [:function, :constant, :type]
```

## Shared solver base

Common infrastructure for the composite solvers: attribute forwarding, sub-solver management, and the
abstraction-based composite template (compute the abstraction, run a control sub-solver, concretize
the controller — each family supplies only the hooks).

```@autodocs
Modules = [Dionysos.Optim]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Continuous-system abstraction solvers

The `Optim.Abstraction` namespace groups the abstraction-based solver families.

```@autodocs
Modules = [Dionysos.Optim.Abstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Uniform grid abstraction

SCOTS-style abstraction on a uniform grid (`GROWTH` / `LINEARIZED` modes), with one control optimizer
per supported specification.

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

Abstraction of timed hybrid systems (per-mode spatial + time abstractions flattened into one
automaton).

```@autodocs
Modules = [Dionysos.Optim.Abstraction.HybridSystemAbstraction]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### PCLF bisimulation quotient

Bisimulation-quotient synthesis for switched systems using a path-complete Lyapunov function, plus
co-safe LTL control on the resulting quotient.

```@autodocs
Modules = [Dionysos.Optim.Abstraction.PCLFBisimulationQuotient]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Trajectory generators

Generators produce a candidate trajectory (open-loop input sequence) for a concrete problem; they are
the building blocks of the lazy solvers and the certifiers.

### Optimizer-based generator

```@autodocs
Modules = [Dionysos.Optim.Abstraction.OptimizerTrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### MPPI generator

```@autodocs
Modules = [Dionysos.Optim.Abstraction.MPPITrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Composite generator

```@autodocs
Modules = [Dionysos.Optim.Abstraction.CompositeTrajectoryGenerator]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Trajectory certifiers

Certifiers take a candidate trajectory and attempt to build a formally-certified tube around it.

### Uniform-grid certifier

```@autodocs
Modules = [Dionysos.Optim.Abstraction.UniformGridTrajectoryCertifier]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Ellipsoidal backward certifier

```@autodocs
Modules = [Dionysos.Optim.Abstraction.EllipsoidalBackwardTrajectoryCertifier]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

### Trajectory certification optimizer

The MOI optimizer that drives a generator into a certifier end to end.

```@autodocs
Modules = [Dionysos.Optim.Abstraction.TrajectoryCertificationOptimizer]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## Discrete-system solvers

Controller synthesis that operates directly on a finite automaton (no abstraction build): worst-case
and optimal cost-to-go fixed-point solvers.

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
