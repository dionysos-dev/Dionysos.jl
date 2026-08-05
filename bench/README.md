# Benchmarks

Performance benchmarks for Dionysos. This is a self-contained environment
([`Project.toml`](Project.toml) sources `Dionysos` via a relative path), so run everything with
`--project=bench`. First time only:

```
julia --project=bench -e 'using Pkg; Pkg.instantiate()'
```

## Cross-solver suite — [`solver_suite.jl`](solver_suite.jl)

Times the full pipeline (setup → abstraction build → controller synthesis) for each
(solver family, problem) pair and writes a comparison table + `results/solver_suite.csv`.

```
julia --project=bench bench/solver_suite.jl [trials] [filter]
```

- `trials` — runs per case, fastest reported (default `1`; use `≥2` to discard compilation latency).
- `filter` — substring of the `"method / problem"` label, e.g. `Gol` or `UniformGrid` (default: all).

The cases live in [`definitions.jl`](definitions.jl); each mirrors the corresponding test under
`test/optim/**`, so they track the live API. Covered families: `UniformGridAbstraction`,
`BemporadMorari`, `BranchAndBound`, `UniformEllipsoidAbstraction`, `LazyEllipsoidsAbstraction`.
Every case is isolated by `try`/`catch`, so one failure is recorded (in the CSV `error` column)
rather than aborting the run. Only open solvers are used (OSQP/HiGHS/Ipopt/Pavito, Clarabel) — no
license-gated ones — so it runs anywhere. Add a case by appending to `CASES`.

## Abstraction-build threading — [`abstraction_threading.jl`](abstraction_threading.jl)

Measures the speed-up of the threaded grid-abstraction backend
(`Symbolic.ThreadedBackend`) over `Symbolic.SequentialBackend` for
`compute_abstract_system_from_concrete_system!`.

```
# one process with N threads:
julia -t <N> --project=bench bench/abstraction_threading.jl [trials] [resolutions] [dt] [du]
# orchestrate several thread counts from one launch:
julia --project=bench bench/abstraction_threading.jl [trials] [resolutions] [dt] [du] [threads]
```

e.g. `julia -t 8 --project=bench bench/abstraction_threading.jl 5 25` runs 5 trials on a 25³ grid.
