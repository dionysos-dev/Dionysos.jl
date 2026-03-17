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
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.build_observation_partition
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.bisimulation_pclf
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.build_sublevel_sequence
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.build_slice_sequence
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.initialize_partitions!
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.initialize_terminal_transitions!
```

```@docs
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient._support_abs_row_on_hyperrectangle
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.gamma_cover_set
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.gamma_cover_region_all_nodes
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.gamma_cover_terminal_all_nodes
Dionysos.Optim.Abstraction.PCLFBisimulationQuotient.build_levels_from_problem
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
