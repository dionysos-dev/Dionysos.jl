# Optim 

This folder contains all the different (abstraction-based or not) solvers that can be used. Note that all the solvers are defined using the MathOptInterface framework as a subtype of  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) by implementig the [`optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.

## Shared solver base

```@docs
Dionysos.Optim.AbstractDionysosOptimizer
Dionysos.Optim.set_field_attribute!
Dionysos.Optim.get_field_attribute
Dionysos.Optim.CompositeDionysosOptimizer
Dionysos.Optim.sub_solvers
Dionysos.Optim.ensure_sub_solvers!
Dionysos.Optim.set_concrete_problem!
```

### Abstraction-based composite template

Shared orchestration for the abstraction-based composite solvers (compute the abstraction, run a
control sub-solver, concretize the controller). Each family supplies only the three hooks.

```@docs
Dionysos.Optim.AbstractionControlOptimizer
Dionysos.Optim.default_abstraction_solver
Dionysos.Optim.build_concrete_controller
Dionysos.Optim.configure_control_solver!
Dionysos.Optim.is_abstraction_computed
```

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
Dionysos.Optim.Abstraction.UniformGridAbstraction.control_solver_for
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerAlternatingSimulationProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerOptimalControlProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerSafetyProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerReachAndStayProblem
Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerCoSafeLTLProblem
```

### Uniform ellipsoid abstraction solver

```@docs
Dionysos.Optim.Abstraction.UniformEllipsoidAbstraction.Optimizer
Dionysos.Optim.Abstraction.UniformEllipsoidAbstraction.control_solver_for
```

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

```@docs
Dionysos.Optim.Abstraction.HybridSystemAbstraction.Optimizer
Dionysos.Optim.Abstraction.HybridSystemAbstraction.control_solver_for
```

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
