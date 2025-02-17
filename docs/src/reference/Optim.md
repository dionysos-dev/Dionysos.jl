# Optim 

This folder contains all the different (abstraction-based or not) solvers that can be used. Note that all the solvers are defined using the MathOptInterface framework: for each solver, we define a subtype of  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) and implement the [`Optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.

## Abstraction-based solvers
```@docs
Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer
Dionysos.Optim.Abstraction.EllipsoidsAbstraction.Optimizer
Dionysos.Optim.Abstraction.HierarchicalAbstraction.Optimizer
Dionysos.Optim.Abstraction.LazyAbstraction.Optimizer
Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer
```

## Other solvers
```@docs
Dionysos.Optim.BemporadMorari.Optimizer
Dionysos.Optim.BranchAndBound.Optimizer
```

## Solvers details

### Uniform Grid Abstraction

```@docs
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_largest_controllable_set
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_largest_invariant_set
```

