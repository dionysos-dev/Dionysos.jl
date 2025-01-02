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
Dionysos.Optim.Abstraction.UniformGridAbstraction._discretize_continuous_system :: Tuple{MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem, Any, Any}
MathOptInterface.optimize! :: Tuple{Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer}
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_continuous_model :: Tuple{Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer}
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_concrete_problem :: Tuple{Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction.build_abstraction :: Tuple{Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_controller_reach! :: Tuple{Any, Any, Any, Vector{Int64}}
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_num_targets_unreachable :: Tuple{Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction._discrete_system :: Tuple{MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem, Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction.compute_controller_safe! :: NTuple{4, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction._maybe_discretized_system :: Tuple{Any, Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_abstract_problem :: Tuple{Dionysos.Problem.SafetyProblem}
Dionysos.Optim.Abstraction.UniformGridAbstraction.solve_abstract_problem :: Tuple{Dionysos.Problem.OptimalControlProblem}
Dionysos.Optim.Abstraction.UniformGridAbstraction._data :: NTuple{4, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_discrete_model :: Tuple{Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer}
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_controller_reach! :: NTuple{7, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction.build_abstract_problem :: Tuple{Any, Dionysos.Symbolic.SymbolicModelList}
Dionysos.Optim.Abstraction.UniformGridAbstraction._validate_model :: Tuple{Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer, Vector{Symbol}}
Dionysos.Optim.Abstraction.UniformGridAbstraction._corresponding_abstract_points :: Tuple{Any, Any, Any}
Dionysos.Optim.Abstraction.UniformGridAbstraction._compute_pairstable :: Tuple{Any, Any}
```

