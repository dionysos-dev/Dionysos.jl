# Problem Types

This module defines a set of structures used to represent different control problems.

All problems are subtypes of the abstract type [`ProblemType`](@ref Dionysos.Problem.ProblemType), which provides a common interface for control problems.

## Supported Problems

- [`EmptyProblem`](@ref Dionysos.Problem.EmptyProblem):  
    Used to construct an abstraction of a dynamical system without solving a control problem.

- [`OptimalControlProblem`](@ref Dionysos.Problem.OptimalControlProblem):  
    A reach-avoid optimal control problem defined over a finite time horizon, supporting state and transition costs.

- [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem):  
    A safety specification problem requiring the system to remain within a safe set for the entire time horizon.

Each of these problem types is detailed below:

```@docs
Dionysos.Problem.ProblemType
Dionysos.Problem.EmptyProblem
Dionysos.Problem.OptimalControlProblem
Dionysos.Problem.SafetyProblem
```