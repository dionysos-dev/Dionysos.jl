# Wrapper

The **JuMP/MOI front-end**: `Model(Dionysos.Optimizer)`. It compiles a JuMP model into a
`MathematicalSystems`/`HybridSystems` object plus a
[`ProblemType`](@ref Dionysos.Problem.ProblemType), then drives an existing solver with them —
it owns no control semantics of its own.

The pipeline is one-directional, `JuMP model → ModelIR → (system, problem) → solver`, with
[`ModelIR`](@ref Dionysos.Wrapper.ModelIR) as the seam: a plain, dependency-free description of
what the user wrote. That is why this module lives in `src/` rather than in an extension —
parsing needs no Symbolics, only compiling a dynamics expression into a callable does (see
[`AbstractDynamicsBackend`](@ref Dionysos.Wrapper.AbstractDynamicsBackend)).

What a model is made of:

- **variables**, whose [`VariableRole`](@ref Dionysos.Wrapper.VariableRole) is inferred from how
  they are used, or declared with [`set_role!`](@ref Dionysos.Wrapper.set_role!);
- **dynamics**, written with `∂`/`Δ` or supplied as a Julia function;
- **specifications**, the markers [`Start`](@ref Dionysos.Wrapper.Start),
  [`Final`](@ref Dionysos.Wrapper.Final), [`Always`](@ref Dionysos.Wrapper.Always) and
  [`EventuallyAlways`](@ref Dionysos.Wrapper.EventuallyAlways), from which the problem type
  follows;
- **modes and transitions** for hybrid models, declared with
  [`@mode`](@ref Dionysos.Wrapper.@mode) and
  [`add_transition!`](@ref Dionysos.Wrapper.add_transition!).

The user guide is `src/wrapper/README.md`; the design and its rationale are in `plan.md`.

## API reference

```@autodocs
Modules = [Dionysos.Wrapper]
Filter  = _is_public
Order   = [:module, :type, :function, :constant, :macro]
```
