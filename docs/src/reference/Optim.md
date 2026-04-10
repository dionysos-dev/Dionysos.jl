# Optim 

This folder contains all the different (abstraction-based or not) solvers that can be used. Note that all the solvers are defined using the MathOptInterface framework as a subtype of  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) by implementig the [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.

## Discrete system solvers
```@docs
Dionysos.Optim.DiscreteSystems.compute_worst_case_cost_controller
Dionysos.Optim.DiscreteSystems.compute_optimal_controller
Dionysos.Optim.DiscreteSystems.compute_worst_case_uniform_cost_controller
```

## Continuous system solvers
### Uniform grid abstraction solver
```@docs
Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerAlternatingSimulationProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerOptimalControlProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerSafetyProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerCoSafeLTLProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.make_out_of_domain_handler
```

### Uniform ellipsoid abstraction solver

### PCLF Bisimulation Quotient solver

```@docs
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.gamma_cover_set
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient._support_abs_row_on_hyperrectangle
```

### Lazy Ellipsoid solver
```@docs
Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer
```

## Hybrid system solvers

## Other solvers

### Bemborad-Morari MIP reformulation

```@docs
Dionysos.Optim.BemporadMorari.Optimizer
Dionysos.Optim.BemporadMorari.julia_function_to_moi
```

### Branch & Bound

```@docs
Dionysos.Optim.BranchAndBound.Optimizer
```
