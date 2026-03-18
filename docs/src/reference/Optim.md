# Optim 

This folder contains all the different (abstraction-based or not) solvers that can be used. Note that all the solvers are defined using the MathOptInterface framework as a subtype of  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) by implementig the [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.

## Abstraction-based solvers
### Uniform grid abstraction solver
```@docs
Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerEmptyProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerOptimalControlProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerSafetyProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerCoSafeLTLProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.make_out_of_domain_handler
```

### Uniform ellipsoid abstraction solver

### Hybrid system abstraction solver

### PCLF Bisimulation Quotient solver

```@docs
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.gamma_cover_set
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient._support_abs_row_on_hyperrectangle
```

### Other abstraction-based solvers
```@docs
Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer
```

## Other solvers
```@docs
Dionysos.Optim.BemporadMorari.Optimizer
Dionysos.Optim.BranchAndBound.Optimizer
```
