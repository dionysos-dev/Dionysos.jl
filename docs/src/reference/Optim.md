# Optim 

This folder contains all the different (abstraction-based or not) solvers that can be used. Note that all the solvers are defined using the MathOptInterface framework as a subtype of  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) by implementig the [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.

## Abstraction-based solvers
### Uniform grid abstraction solver
```@docs
Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerEmptyProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerOptimalControlProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerSafetyProblem
```

### Other abstraction-based solvers
```@docs
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

## Timed Hybrid Abstraction
```@docs
Dionysos.Optim.Abstraction.TimedHybridAbstraction.TimedHybridProblemSpecs
Dionysos.Optim.Abstraction.TimedHybridAbstraction.TimedHybridOptimalControlProblem
Dionysos.Optim.Abstraction.TimedHybridAbstraction.TimedHybridSafetyProblem
Dionysos.Optim.Abstraction.TimedHybridAbstraction.solve_timed_hybrid_problem
Dionysos.Optim.Abstraction.TimedHybridAbstraction.get_closed_loop_trajectory
Dionysos.Optim.Abstraction.TimedHybridAbstraction.Optimizer
```
