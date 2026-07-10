# Coding conventions

This page is the **single source of truth** for how Dionysos code is written. It exists so that every
module looks and behaves the same way, so that refactors land against a fixed target, and so that
contributors (and AI assistants) don't have to reverse-engineer the house style from the surrounding
code. New code must follow it; code being touched should be migrated toward it.

These conventions are enforced softly by review and by [`Aqua`](https://github.com/JuliaTesting/Aqua.jl)
+ [`JuliaFormatter`](https://github.com/domluna/JuliaFormatter.jl) (see [Set up](@ref) and the
`.JuliaFormatter.toml` at the repo root). Format every file before committing.

## 1. Interfaces and required methods

Abstract types declare a *contract*: a set of methods each concrete subtype must implement. Declare each
required method **once**, and make an unimplemented method **fail loudly**. Never leave a silent
empty-body stub — a method that returns `nothing` when unimplemented turns a missing-method bug into a
silent wrong-answer bug.

```julia
abstract type AbstractWidget end

# Required — every subtype must implement these:
get_size(w::AbstractWidget) = error("get_size not implemented for $(typeof(w))")
compute!(w::AbstractWidget) = error("compute! not implemented for $(typeof(w))")

# Derived — implemented once in terms of the required methods:
is_empty(w::AbstractWidget) = get_size(w) == 0
```

Rules:
- Use the message form `error("<name> not implemented for $(typeof(x))")`. Do **not** use empty bodies
  (`function f(::AbstractWidget) end`), and standardize the wording (not `"implement f"` /
  `"Not implemented"` / `"not implemented"` mixed together).
- Document, next to the abstract type, which methods are **required** vs **derived**.
- If a subtype extends an interface owned by another package (e.g. `HybridSystems.add_state!`), extend
  *that* function — do not declare a second same-named function in a Dionysos module, or the two will not
  dispatch to each other.

## 2. Type stability

Type stability is a hard requirement, especially on the abstraction / automaton / nearest-neighbour hot
paths.

- **No `::Any` struct fields.** Parameterize the struct instead.
- **No bare `::Function` struct fields.** Make the callable a type parameter (`struct S{F}; f::F; end`);
  reach for `FunctionWrappers.jl` only when a concrete callable type genuinely cannot be a parameter.
- **Parametric structs must bind every field type.** `Tree{S}` with `root::Node{S}`, never `Tree` with
  `root::Node` (an unbound `UnionAll` field forces dynamic dispatch on every access).
- Avoid `[]` (`Vector{Any}`) accumulators; write `T[]` with a concrete `T`.
- Prefer `Union{Nothing, T}` only for genuinely optional state; a field that is "always present after
  construction" should not be nullable.
- Check hot paths with `@code_warntype` / [`JET.jl`](https://github.com/aviatesk/JET.jl) and cover them
  with `@inferred` tests.

## 3. Naming

- **Modules**: `CamelCase`. **Concrete types**: `CamelCase`. **Abstract types**: `Abstract`-prefixed
  (`AbstractMapping`, `AbstractController`).
- **Functions**: `snake_case`. No `camelCase` — rename on sight (`centerDistance → center_distance`,
  `kNearestNeighbors → k_nearest_neighbors`, `get_nNodes → get_n_nodes`, `get_max_Node → get_max_node`).
- Accessors: `get_<noun>` / `get_<noun>_by_<key>`; enumerators: `enum_<plural>`; predicates:
  `is_<adj>` / `has_<noun>`. Counts: `num_<plural>` (or `get_n_<noun>` — pick one per subsystem and be
  consistent).
- **Mutating functions end with `!`**; private helpers start with `_`.
- **One accessor per concept.** Do not ship both `get_dim` and `get_dims`, or both `volume` and
  `get_volume`, for the same quantity. If two functions have identical bodies, keep one.
- Unicode identifiers are idiomatic where they match the maths (`∂`, `Δ`, `α`, `δx`); `∈`/`∉` for
  membership.

## 4. Verbosity and logging

- **No `print` / `println` in library code.** Use `@info` / `@warn` / `@debug`, or
  [`ProgressMeter`](https://github.com/timholy/ProgressMeter.jl) (already a dependency) for progress
  bars. Gate repeated warnings with `maxlog`.
- **One verbosity knob**: an integer `print_level` (`0` = silent, `1` = default, `2` = detailed) plus
  support for `MOI.Silent`. Do not introduce `log_level`, `verbose`, `debug`, etc. in new code.

## 5. Modules, imports, and comments

- **`import`, don't `using`.** Prefer `import Module` (or `import ..Utils as UT`) and call through the
  prefix — `LazySets.dim(x)`, `UT.set_union(v)` — so every name's origin is explicit. Reserve `using`
  for the rare need to pull in operators/macros that must be unqualified (e.g. `using LinearAlgebra` for
  `I` and `\`), and keep it to that.
- **Comment the *why*, not the *what*.** Do not add comments that restate the code; keep only
  non-obvious rationale, invariants, gotchas, or correctness constraints. Docstrings cover the API.

## 6. Immutability

- Value objects (problems, sets, approximations, trajectories) are **immutable `struct`s**. Build a new
  object rather than mutating in place.
- Use `mutable struct` only where in-place mutation is part of the design (incremental builders, indexes,
  caches, solver state). State the reason in a comment when it is not obvious.

## 7. Reuse the ecosystem (behind an interface)

Prefer a mature package over a hand-rolled equivalent, wrapped behind a small Dionysos interface so the
implementation can be swapped:

| Need | Use |
| :-- | :-- |
| Sets & set algebra (∪, ∩, `\`, membership, boxes, ellipsoids, polytopes) | `LazySets` |
| Sorted containers, heaps, priority queues | `DataStructures` |
| Nearest-neighbour queries | `NearestNeighbors` |
| Small labelled graphs | `Graphs` |
| Special functions (e.g. Γ for n-ball volume) | `SpecialFunctions` |
| System / hybrid-system representations | `MathematicalSystems` / `HybridSystems` |

Add a dependency **in the same change that first uses it** — `Aqua`'s stale-dependency check fails if a
package is declared in `[deps]` but never `import`ed/`using`ed.

## 8. Solvers (`Optim`)

Every solver is a `MOI.AbstractOptimizer`. To keep the "swap / compare / benchmark the solver" story
first-class:

- Build on the shared optimizer base rather than re-implementing the `MOI.set` / `MOI.get` /
  `MOI.is_empty` / `reset!` plumbing by hand.
- Prefer typed options (or validated `RawOptimizerAttribute`s) over unchecked string keys; a wrong key or
  wrong value type should error at set time, not deep inside `optimize!`.
- Use the standard field names `print_level` and `solve_time_sec`.
- Select behaviour by multiple dispatch on the problem/system type, not by an `isa`-chain edited in
  several places.

## 9. Tests, formatting, docs

- Run the suite: `julia --project -e 'using Pkg; Pkg.test()'`.
- Format before committing: `julia -e 'using JuliaFormatter; format(".")'` (CI fails on any diff).
- Every **exported** symbol needs a docstring (the docs build errors otherwise).
- Behaviour-preserving refactors must keep the golden-output regression suite (`test/regression/`, once
  it exists) unchanged.
