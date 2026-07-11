# CLAUDE.md

Operating guide for AI coding agents (Claude Code) working in **Dionysos.jl**. Read this before
making changes: it explains what the project is, how the code is organized, the one architectural
contract that matters most (the solver interface), the Julia conventions you must imitate, the exact
commands to run, and the repo-specific traps to avoid.

---

## 1. What Dionysos is & why it exists

Dionysos is a **framework for symbolic (abstraction-based) control**. It is the software of the ERC
project *Learning to Control* (L2C).

**The mission.** Today, controlling a complex system usually means a team of expert engineers, each
with deep knowledge of the plant, hand-crafting an ad hoc controller over months at significant cost.
Dionysos aims to change that paradigm into an automatic pipeline:

> **describe the system → select the problem (specification) → pick a solver → Dionysos automatically
> computes a controller with a certificate (formal guarantee).**

The goal is to make correct-by-construction control synthesis *accessible* — including to small
companies that have no dedicated controls/IT team — and to drastically cut controller-design time.

**Symbolic control 101.** The system is *abstracted* into a finite-state automaton (a **symbolic
model**) by discretizing its variables. Working on this finite object lets a computer synthesize
controllers systematically, even for complex specifications (e.g. LTL), using graph algorithms
(Dijkstra, A\*, fixed-point iterations). The price:

- **Curse of dimensionality** — the number of abstract states grows exponentially with the state
  dimension.
- **Non-determinism from discretization** — over-approximating a cell's dynamics introduces spurious
  transitions, which can make the abstract problem *infeasible* even when the concrete one is not.

The mitigation, and a core research direction of the toolbox, is **smart / lazy abstractions**: instead
of an a-priori uniform grid, co-design the abstraction and the controller and compute only the part of
the abstraction that is actually needed.

**Framing for an agent.** Dionysos is an **ecosystem, not one algorithm**. Its value is a *common
interface* that lets you swap the solver and then test / compare / benchmark algorithms on the same
control problem. Preserving that interface is the single most important thing when extending the code.

Background reading: [docs/src/manual/manual.md](docs/src/manual/manual.md) and
[docs/src/manual/abstraction-based-control.md](docs/src/manual/abstraction-based-control.md).

---

## 2. The framework in one paragraph

A control problem is a pair `(𝒮, Σ)`: a **system** `𝒮` (a
[`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl) or
[`HybridSystems`](https://github.com/blegat/HybridSystems.jl) object) plus a **specification** `Σ` (a
`Dionysos.Problem.ProblemType`). It is solved by a **solver** `𝒪` that is a
`MathOptInterface.AbstractOptimizer`, driven through the JuMP/MOI interface. This MOI contract is the
architectural keystone: **every algorithm is a swappable `Optimizer`**, so a task can be re-solved,
compared, and benchmarked by swapping the optimizer rather than rewriting the model. Abstraction-based
solvers turn the infinite-state system into a finite automaton via discretization, synthesize a
correct-by-construction controller on it, then *concretize* it back to the original system.

---

## 3. Repository map (root level)

| Path | What it is |
| :--- | :--- |
| [`src/`](src/) | The `Dionysos` package (all library code). |
| [`ext/`](ext/) | Package extensions — optional-dependency glue (Plots, Symbolics, Spot, CSV, RigidBodyDynamics). |
| [`test/`](test/) | Test suite; mirrors `src/` layout. Entry point [`test/runtests.jl`](test/runtests.jl). |
| [`docs/`](docs/) | Documenter.jl site + Literate.jl examples. Build script [`docs/make.jl`](docs/make.jl). |
| [`problems/`](problems/) | Reusable **benchmark problem library** (e.g. path planning, DC-DC, pendulum). |
| `utils/` (root) | **Runnable case-study scripts**, grouped by example — **not** library code. |
| [`bench/`](bench/) | Benchmarks (BenchmarkTools). |
| `control_server/`, `BipedRobot/`, `paper/`, `assets/` | Auxiliary app, robot demo, paper artifacts, images. |

> ⚠️ The root `utils/` folder (example drivers) is **not** the same as the `src/utils/` module (the
> `Utils` library). Don't confuse them.

---

## 4. Core architecture — the six modules

Top-level module [`src/Dionysos.jl`](src/Dionysos.jl) includes six submodules **in dependency order**:

```
Utils → System → Problem → Mapping → Symbolic → Optim
```

Reuse these **standard aliases** everywhere (they are established at the top of each module — match
them exactly):

```julia
const UT = Utils    # Dionysos.Utils
const ST = System   # Dionysos.System
const PR = Problem  # Dionysos.Problem
const MP = Mapping  # Dionysos.Mapping
const SY = Symbolic # Dionysos.Symbolic
const OP = Optim    # Dionysos.Optim
const AB = OP.Abstraction
```

| Module | File | Responsibility |
| :--- | :--- | :--- |
| `Utils` | [src/utils/utils.jl](src/utils/utils.jl) | Foundational helpers on top of LazySets: sets *are* LazySets (`UT.box`/`LazySets.Hyperrectangle`, `LazySets.Ellipsoid` with shape matrix `Q = P⁻¹`, unions/set-minus), plus what LazySets lacks — exact ellipsoid kernels (`is_included`/`is_disjoint`), periodic splitting (`set_in_period`), quadratic-form bridge (`get_quadratic_form`) — callable cost functions (`functions.jl`), data structures, `search/RRT.jl`, scalar optimization (`numeric/`), plotting recipes, PCLF. |
| `System` | [src/system/system.jl](src/system/system.jl) | Concrete dynamical systems (extends MathematicalSystems / HybridSystems), time discretization, reachable-set over/under-approximations of the dynamics (`approximation/`, with the `input_cache`/`reach_set` per-input hoisting hooks the symbolic backend builds on), local affine approximations (`build_affine_approximation` + providers), controllers (trait-based protocol, plain-data/serializable), trajectories + the closed-loop simulation engine, ellipsoidal transition synthesis (`solve_transition`). |
| `Problem` | [src/problem/problems.jl](src/problem/problems.jl) | Control-task **specifications** = `ProblemType` subtypes (solver-independent). |
| `Mapping` | [src/mapping/mapping.jl](src/mapping/mapping.jl) | Concrete ↔ abstract **discretization**: grids, cells, `AbstractMapping`, inclusion modes `INNER/OUTER/CENTER`. |
| `Symbolic` | [src/symbolic/symbolic_model.jl](src/symbolic/symbolic_model.jl) | Builds the finite **automaton abstraction** (`SymbolicModel`) from a system + mapping; parallel build backends (threaded / distributed / SLURM). |
| `Optim` | [src/optim/optim.jl](src/optim/optim.jl) | The **solver catalog**. |

Also in `src/Dionysos.jl`: the JuMP `NonlinearOperator`s `∂`, `Δ`, `final`, `start` (used to express
dynamics and reach/target constraints in a JuMP model), plus stub functions for extension-only features.

---

## 5. The solver contract (the keystone — read this before touching `src/optim/`)

**Every solver is a submodule exposing an `Optimizer <: MOI.AbstractOptimizer`.** It is:

1. Instantiated: `optimizer = MOI.instantiate(SomeFamily.Optimizer)`
2. Configured: `MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)` (and
   `"state_grid"`, `"input_grid"`, `"time_step"`, `"print_level"`, …)
3. Run: `MOI.optimize!(optimizer)`
4. Queried: `MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))`

**Conceptual pipeline** (canonical: [uniform_grid_abstraction.jl](src/optim/continuous_systems/UniformGridAbstraction/uniform_grid_abstraction.jl)):

```
concrete problem → abstraction (symbolic model) → abstract problem
                 → abstract controller → concrete controller → closed-loop trajectory
```

The abstraction is **cached**: switching the control task (e.g. safety → reachability) on the same
system does *not* recompute it. Solvers **compose** — a high-level optimizer holds sub-solvers (an
`abstraction_solver` and a `control_solver`) and forwards attribute `set`/`get` to them.

### Two user entry styles

**(a) JuMP model — canonical**, see [docs/src/examples/solvers/Path planning.jl](docs/src/examples/solvers/Path%20planning.jl):

```julia
using Dionysos, JuMP, StaticArrays
model = Model(Dionysos.Optimizer)
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = x_start[i])
@variable(model, -1 <= u[1:2] <= 1)
@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))   # dynamics via ∂
@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))       # target set
set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", Dionysos.Mapping.GridFree(x0, hx))
optimize!(model)
concrete_controller = get_attribute(model, "concrete_controller")
```

**(b) Direct MOI** on a specific family optimizer, see [docs/src/examples/solvers/Lazy Ellipsoids Abstraction.jl](docs/src/examples/solvers/Lazy%20Ellipsoids%20Abstraction.jl).

### Implemented solver families

| Family | Location |
| :--- | :--- |
| Uniform grid abstraction (SCOTS-style; `GROWTH`/`LINEARIZED`) | [src/optim/continuous_systems/UniformGridAbstraction/](src/optim/continuous_systems/UniformGridAbstraction/) |
| Uniform ellipsoid abstraction | [src/optim/continuous_systems/UniformEllipsoidAbstraction/](src/optim/continuous_systems/UniformEllipsoidAbstraction/) |
| Lazy ellipsoids abstraction (RRT + SDP/Lyapunov) | [src/optim/continuous_systems/lazy_ellipsoids_abstraction.jl](src/optim/continuous_systems/lazy_ellipsoids_abstraction.jl) |
| Hybrid-system abstraction | [src/optim/hybrid_systems/HybridSystemAbstraction/](src/optim/hybrid_systems/HybridSystemAbstraction/) |
| PCLF bisimulation quotient | [src/optim/hybrid_systems/PCLFBisimulationQuotient/](src/optim/hybrid_systems/PCLFBisimulationQuotient/) |
| Discrete-system solvers (operate directly on an automaton) | [src/optim/discrete_systems/](src/optim/discrete_systems/) |
| Trajectory generators (MPPI, optimizer-based, composite) | [src/optim/trajectory_generators/](src/optim/trajectory_generators/) |
| Trajectory certifiers | [src/optim/trajectory_certifiers/](src/optim/trajectory_certifiers/) |
| Bemporad–Morari (MIQP for PWA) | [src/optim/bemporad_morari.jl](src/optim/bemporad_morari.jl) |
| Branch and bound | [src/optim/branch_and_bound.jl](src/optim/branch_and_bound.jl) |
| Q-learning | [src/optim/q_learning.jl](src/optim/q_learning.jl) |

### Checklist: adding a new solver

1. Create a submodule under `src/optim/…` exposing `Optimizer <: MOI.AbstractOptimizer`.
2. Implement `MOI.is_empty`, `MOI.set`, `MOI.get` (raw attributes), and `MOI.optimize!`.
3. Wire the `include(...)` into [src/optim/optim.jl](src/optim/optim.jl).
4. Add a **docstring** for every exported symbol — the docs build errors on missing docstrings.
5. Add a test under `test/optim/…` and wire it into [test/runtests.jl](test/runtests.jl).
6. If user-facing, add a Literate example under `docs/src/examples/solvers/`.

---

## 6. Domain vocabulary (glossary)

- **Concrete vs abstract** — everything comes in pairs: `concrete_problem`/`abstract_problem`,
  `concrete_controller`/`abstract_controller`, `concrete_system`/`abstract_system`. Conversion helpers:
  `get_concrete_state`/`get_abstract_state`, `quantize_controller`.
- **Concrete system** — the real system `ẋ = f(x,u)` (a `MathematicalSystems` object). `system.X` is its
  state domain.
- **Symbolic model / automaton** — a finite, *sound* over-approximation of the concrete system,
  `SymbolicModel{N,M}` (concretely `SymbolicModelList`) wrapping an automaton with transitions
  `(target, source, symbol)`. `N` = state dim, `M` = input dim. See
  [src/symbolic/symbolic_model.jl](src/symbolic/symbolic_model.jl).
- **Grid / cell / mapping / domain** — the space is discretized by a **grid** `GridFree(x0, hx)` (center
  `x0`, step `hx`); a **cell** is an integer tuple `pos` with coordinate `x0 + hx .* pos`; a **Mapping**
  (`AbstractMapping`, e.g. `ExplicitGridMapping`) is the bijection between integer state labels `1:n` and
  cells; a `MappingSet` is a collection of cells with inclusion mode `INNER`/`OUTER`/`CENTER`.
- **Specification / problem** — a `ProblemType` subtype: `OptimalControlProblem` (reach-avoid),
  `SafetyProblem`, `ReachAndStayProblem`, `CoSafeLTLProblem`, and the abstraction-only
  `AlternatingSimulationProblem` / `BisimulationQuotientProblem`. See
  [src/problem/problems.jl](src/problem/problems.jl). Infinite horizons use the `Infinity()` sentinel.
- **Growth bound** — a function over-approximating a cell's reachable set (to build sound transitions);
  supplied via `ST.ContinuousTimeGrowthBound(system; jacobian_bound)`. `approx_mode`
  selects `GROWTH` vs `LINEARIZED`.
- **Post / pre** — `post(sym, q, u)` = successor abstract states; `pre(sym, target)` = predecessors.
  Fixed-point pre-image computation drives reachability synthesis.
- **Alternating simulation / bisimulation** — the soundness relations guaranteeing a controller
  synthesized on the abstraction is valid on the concrete system.
- **Value / Lyapunov function** — the certificate: `abstract_value_function` (worst-case cost-to-go) and
  the Lyapunov functions for the ellipsoid abstractions.
- **Controller / feedback** — `AbstractController` splits into **static** (state-feedback map `q ↦ u`)
  and **dynamic** (has internal memory). Protocol: `initial_state`, `update_state`, `output_control`,
  `is_defined`, `domain`. See [src/system/controllers/controller.jl](src/system/controllers/controller.jl).

---

## 7. Julia conventions & idioms (imitate these)

> Full, authoritative version: [docs/src/developers/conventions.md](docs/src/developers/conventions.md)
> — the house style for interfaces, type stability, naming, logging, immutability, ecosystem reuse, and
> solver structure. The summary below is the short form.

**Naming**
- Modules: `CamelCase`, one per subsystem. Abstract types: `Abstract`-prefixed
  (`AbstractMapping`, `AbstractController`). Concrete types: `CamelCase`.
- Functions: `snake_case`. Getters `get_<noun>_by_<key>`, enumerators `enum_<plural>`, predicates
  `is_<adj>`/`has_<noun>`. **Mutating functions end with `!`** (`add_transition!`, `compute_post!`).
- Private helpers: leading underscore (`_invInclMode`, `_diff`).
- **Unicode identifiers are idiomatic**: `∂`, `Δ`, `α`, `δx`, `q′`, and `∈`/`∉` for set membership.

**Interface pattern** — the dominant idiom. Define an abstract type plus method *stubs* that
`error("implement …")` or have empty bodies; concrete subtypes override them. Behavior is selected by
**multiple dispatch** (a generic fallback + one method per type), e.g. `discretize_problem` /
`trajectory_success` in `problems.jl`.

**Types & data**
- **Parametric structs carry dimensions and field types** for type stability:
  `SymbolicModel{N,M}`, `AbstractMapping{N,T}`, `OptimalControlProblem{S,XI,XT,XC,T}`.
- **StaticArrays everywhere** for small fixed-size data — `SVector`/`SMatrix` for coordinates, grid
  steps, system matrices, ellipsoid shapes.
- Validating/normalizing constructors (e.g. `UT.box` rejects crossed bounds; an empty
  region is `LazySets.EmptySet`, never a sentinel).
- Sentinel/singleton types and `@enum` instead of magic values (`struct Infinity end`,
  `@enum INCL_MODE INNER OUTER CENTER`).
- `nothing` + short-circuit guards in hot paths (`q in domain(ctrl) || return false`).
- Plot logic lives in `@recipe`/`@series` functions next to each type — keep it out of core logic.

**Reuse before writing.** Prefer existing primitives — `UT.box` / `LazySets.Hyperrectangle`,
`LazySets.Ellipsoid` (construct from a quadratic-form matrix `P` via `Q = inv(P)`; never swap
the two silently), `MP.GridFree`, the LazySets API (`center`, `low`/`high`, `∈`,
`box_approximation`, `vertices_list`, `sample`) and the data structures in `src/utils/` — over
re-implementing. Set predicates go through `UT.is_included` / `UT.is_disjoint` (Dionysos-owned
verbs; `Base` methods on two LazySets-owned types are piracy and fail the Aqua gate). Grid
discretization (`MP.get_states_from_set`) accepts any bounded `LazySet` — zonotopes, balls,
polytopes included.

---

## 8. Performance guidance

- **Type stability** via parametric fields; keep functions inferrable.
- **StaticArrays** for small vectors/matrices (avoids heap allocation).
- Prefer **generators** (`enum_states`, `enum_coords`) over materializing arrays.
- `@inline` small helpers; `sizehint!` before push loops; `get!(dict, key) do … end` for memoized id
  assignment.
- `@inbounds`/`@simd` are used **sparingly** in verified hot loops (automaton/grid), *not* pervasively —
  don't scatter them.
- For large abstractions, use the parallel build backends in
  [src/symbolic/](src/symbolic/) (threaded / distributed / SLURM) rather than the sequential path.

---

## 9. Developer workflow & commands

Julia ≥ 1.10 (CI tests 1.10 and current `1`). Commands assume the repo root.

**Run tests — don't run the whole suite for every change; it's slow.** Prefer the narrowest scope:

1. **Run only the test file(s) you touched, standalone.** Every test file is a self-contained module
   (its own `using`) and *must be runnable on its own* in the test environment
   ([test/Project.toml](test/Project.toml), which sources `Dionysos` via a relative path):
   ```
   julia --project=test test/optim/UniformGridAbstraction/unit_test_reachability.jl
   ```
   First time only: `julia --project=test -e 'using Pkg; Pkg.instantiate()'`. If a file is *not*
   standalone-runnable (e.g. a missing import), fix it — add the missing `using`/`import`.
2. **Fast subset smoke check:** `julia --project -e 'using Pkg; Pkg.test(; test_args = ["--fast"])'`
   runs everything except suites tagged `:slow` in [test/runtests.jl](test/runtests.jl) (heavy
   end-to-end solver pipelines + Aqua). The driver prints per-file timings and a slowest-first summary.
3. **Full gate before committing / in CI:** `julia --project -e 'using Pkg; Pkg.test()'`.

New test files must be added to the `TEST_FILES` list in [test/runtests.jl](test/runtests.jl) (tag a
slow suite `:slow`). **Coverage check:** to confirm new lines are exercised, put `@show @__LINE__`
around them and look for those line numbers in the test log — if they don't appear, the lines are
untested.

**Persistent Julia REPL (avoid precompilation).** Julia pays a large precompilation cost per process.
For quick iterations (scripts under ~30 s), keep a long-lived REPL in `tmux` and send commands to it
instead of relaunching `julia`:
```
tmux new-session -d -s julia -x 220 -y 50
tmux send-keys -t julia 'julia --project=test' Enter          # then wait for the julia> prompt
tmux send-keys -t julia 'using Pkg; Pkg.instantiate()' Enter  # first time only
tmux send-keys -t julia 'using Revise' Enter                  # hot-reload edited code
tmux send-keys -t julia 'includet("test/…/foo.jl")' Enter     # includet = tracked by Revise
# wait for completion, then inspect:
while ! tmux capture-pane -t julia -p | grep -q "^julia>"; do sleep 1; done
tmux capture-pane -t julia -p
```
Tear down with `tmux kill-session -t julia`. (Unix/tmux; on Windows use a persistent VSCode Julia REPL.)

**Julia & dev packages.** Julia is managed by `juliaup` (`~/.julia/juliaup`); local dev checkouts of
packages live in `~/.julia/dev`. If you need the source of a package that isn't there, ask for it to be
cloned. When `Pkg.develop`-ing a package use a **relative path**, never an absolute one — absolute paths
(`/home/<user>/…`) break the shared `Manifest.toml` across machines.

**Format — REQUIRED before every commit** (CI `format_check.yml` fails on any diff):
```
julia -e 'using JuliaFormatter; format(".")'
```
Rules in [.JuliaFormatter.toml](.JuliaFormatter.toml): `always_use_return = true` (every function ends
with an explicit `return`), `separate_kwargs_with_semicolon`, `always_for_in`, `whitespace_in_kwargs`,
`whitespace_typedefs`, `surround_whereop_typeparameters`, `indent_submodule`, 4-space indent.

**Build docs**
```
julia --project=docs/ docs/make.jl
```
Documenter + Literate + DocumenterCitations. `makedocs(modules=[Dionysos])` **errors on any exported
symbol lacking a docstring** — so new exported symbols must be documented. For fast local iteration,
comment out the Literate loop in [docs/make.jl](docs/make.jl). Do **not** commit `docs/Project.toml`
changes produced by `dev .`.

**CI** ([.github/workflows/](.github/workflows/)): `ci.yml` (test matrix 1.10 + `1`, coverage),
`format_check.yml`, `aqua.yml`, `documentation.yml`, `BipedRobot.yml`.

**Git workflow** ([docs/src/developers/git.md](docs/src/developers/git.md)): never commit to `master`;
branch per change; **format before committing**; open a PR (CI runs tests + format check). The default
branch for PRs is `master`.

**Package extensions.** Optional features live in [ext/](ext/) behind `[weakdeps]`/`[extensions]` in
[Project.toml](Project.toml): `DionysosPlotsExt`, `DionysosSymbolicsExt` (Symbolics + IntervalArithmetic),
`DionysosSpotExt`, `DionysosCSVExt` (CSV + DataFrames), `DionysosMathOptSymbolicADExt`,
`DionysosRigidBodyDynamicsExt`. The stub functions in `src/Dionysos.jl` error until the corresponding
extension is loaded (`using Plots`, `using Symbolics`, …).

---

## 10. Gotchas / repo-specific traps

- **Stale `Domain` naming.** Some docstrings/comments say `Dionysos.Domain.INNER`, but the live module
  is **`Dionysos.Mapping`** (`INNER`/`OUTER`/`CENTER`). "Domain" is a former name — verify against the
  code before trusting a comment.
- **Two `utils/`.** Root `utils/` = case-study scripts; `src/utils/` = the `Utils` library module.
- **Wire new tests in.** [test/runtests.jl](test/runtests.jl) is a `TEST_FILES` list (path + optional
  `:slow` tag) run in a timed loop; its include paths are slightly flatter than the deeply nested `src/`
  tree. Add every new test file to that list, and mirror the *actual* `src/` layout for source.
- **Cross-platform.** The repo is developed both in a Linux container and on Windows. The `julia
  --project …` commands are portable; the `tmux` persistent-REPL recipe is Unix-only (on Windows use a
  persistent VSCode Julia REPL instead).
