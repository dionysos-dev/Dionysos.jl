# Problem 

This folder contains structures that are used to encode which kind of problem you want to solve. Each problem is encoded as a ProblemType.

```@docs
Dionysos.Problem.ProblemType
```

So far, two types of problems have been considered: the reach-avoid optimal control problems and the safety control problems.

```@docs
Dionysos.Problem.EmptyProblem
Dionysos.Problem.OptimalControlProblem
Dionysos.Problem.SafetyProblem
```