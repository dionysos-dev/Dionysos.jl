# Problem

Solver-independent **specifications**. Every specification is a subtype of
[`ProblemType`](@ref Dionysos.Problem.ProblemType), split into
[`ControlProblem`](@ref Dionysos.Problem.ControlProblem) (a controller is synthesized: reach-avoid,
safety, reach-and-stay, co-safe LTL) and
[`AbstractionProblem`](@ref Dionysos.Problem.AbstractionProblem) (no control objective; parametrizes a
reusable abstraction). Infinite horizons use the [`Infinity`](@ref Dionysos.Problem.Infinity)
sentinel.

## API reference

```@autodocs
Modules = [Dionysos.Problem]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
