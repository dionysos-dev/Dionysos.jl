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
Dionysos.Optim.Abstraction.UniformGridAbstraction._get_domain_list
Dionysos.Optim.Abstraction.UniformGridAbstraction._discretize_continuous_system
MathOptInterface.optimize!
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_continuous_model
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_concrete_problem
Dionysos.Optim.Abstraction.UniformGridAbstraction.build_abstraction
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_controller_reach!
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_num_targets_unreachable
Dionysos.Optim.Abstraction.UniformGridAbstraction._discrete_system
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_controller_safe!
Dionysos.Optim.Abstraction.UniformGridAbstraction._maybe_discretized_system
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_abstract_problem
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_abstract_problem
Dionysos.Optim.Abstraction.UniformGridAbstraction._data
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_discrete_model
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_controller_reach!
Dionysos.Optim.Abstraction.UniformGridAbstraction.build_abstract_problem
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_model
Dionysos.Optim.Abstraction.UniformGridAbstraction._corresponding_abstract_points
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_pairstable
```

