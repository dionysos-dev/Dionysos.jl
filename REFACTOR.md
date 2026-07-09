# Dionysos.jl v0.2 refactor roadmap

Living roadmap for the whole-toolbox refactor toward **v0.2**: modularity, consistency, type-stability,
and an extensible solver layer. House style lives in
[docs/src/developers/conventions.md](docs/src/developers/conventions.md); this file is the *plan and
status*.

## Goal

A modular, consistent, efficient, extensible toolbox — clean per-module APIs, ecosystem reuse behind
clean interfaces, type-stable hot paths, and a plugin-style solver layer that makes "change / compare /
benchmark the solver" first-class.

## Confirmed decisions

1. **Ecosystem = interface-first.** Define clean internal interfaces (Set, Automaton, NearestNeighbor,
   Optimizer); back them with the ecosystem where it clearly wins — finish the LazySets migration;
   adopt Graphs + SpecialFunctions (already transitive deps), DataStructures, NearestNeighbors. Keep a
   **purpose-built CSR automaton** behind the interface for the million-edge hot path.
2. **Bold clean-break → v0.2.** Freely redesign/rename/split internal types. Keep the **public surface**
   working as the regression harness: the JuMP `Dionysos.Optimizer` DSL (`∂/Δ/start/final`, `∉`), the
   `problems/*` benchmark library, and `docs/src/examples/*`. Ship v0.2 + `MIGRATION.md`.
3. **Foundations first, then bottom-up.** Phase 0 (conventions + safety net + shared bases + bug fixes),
   then Utils → System → Mapping → Problem → Symbolic → Optim.

## Why (audit findings that motivate this)

- **Reinvented ecosystem code, half-migrated.** A hand-rolled `AbstractLazySet` +
  `LazySetUnion/Intersection/Minus` that collides with `LazySets`' own names; hand-rolled `Tree` +
  brute-force kNN (O(n log n) per RRT query), `SortedTupleSet`, a custom n-ball Γ, and a hand-built
  labeled-multigraph automaton — all duplicating LazySets / NearestNeighbors / DataStructures /
  SpecialFunctions / Graphs.
- **No uniform interface style** — four "required method" conventions coexist, several failing *silently*.
- **Pervasive type-instability** — `SymbolicSystem` (20 `::Any` fields), `Tree` (unbound UnionAll),
  `::Function`/`::Any` fields, `transition_cost::Any`, 12-param `SymbolicModelList`.
- **Solver boilerplate + fragile dispatch** — ~15 identical `MOI.set/get` bodies, 3 drifted ~80-line
  wrappers, no shared base, string-keyed attrs without validation, `isa`-chain dispatch (3 edit sites).
- **Leaky/duplicated modeling** — "domain" means three things; state-sets aren't self-contained; 5+
  trajectory types (one used); static/dynamic controllers are naming not types; problems needlessly
  `mutable`; `Domain → Mapping` rename incomplete.
- **~10 latent bugs** in dead/untested paths (fixed in Phase 0).

## Guiding principles (see conventions.md)

One interface convention (loud `error`, no silent stubs) · type stability is a hard rule (no `::Any`/
bare `::Function` fields) · snake_case + one accessor per concept · one verbosity knob (`print_level`,
no `print`/`println` in library code) · immutable value objects · reuse the ecosystem behind interfaces.

## Phased roadmap & status

### Phase 0 — Foundations & safety net ✅ DONE
- ✅ ~11 latent bug fixes (dead/untested paths).
- ✅ Conventions guide + CLAUDE.md workflow.
- ✅ Timed, taggable test driver (`test/runtests.jl` `TEST_FILES` list) with `--fast` opt-out
  (`Pkg.test(; test_args=["--fast"])` skips `:slow` suites).
- ✅ Golden-output regression harness (`test/regression/`) — pins the flagship UniformGrid
  path-planning fingerprint (states 26285, transitions 9705631, controllable 19703).
- ✅ Shared `AbstractDionysosOptimizer` base (`src/optim/optimizer_common.jl`) — validated field-backed
  `RawOptimizerAttribute` get/set + `SolveTimeSec`; discrete_systems group migrated. (Adoption across
  the remaining ~15 optimizers happens in Phase 6.)
- ⏸ Dep promotion deferred to first use (Aqua stale-deps blocks unused deps).

### Phase 1 — Utils (biggest ecosystem win) ⬜ IN PROGRESS
- Sets → **LazySets** as the canonical representation: delete the hand-rolled lazy-set algebra
  (`src/utils/sets/lazy_set_operations.jl`) and custom set types; route through `UnionSetArray` /
  `Intersection(Array)` / `Complement`. Keep Dionysos-only types (`DeformedRectangle`, integer `IndexBox`
  for grids) behind the same interface. `SpecialFunctions.gamma` for ellipsoid volume. Collapse duplicate
  accessors (`get_dim`/`get_dims`, `volume`/`get_volume`, `expand≡get_sublevel_set≡*`,
  `transform≡affine_transformation`).
- Data structures & search: type-stable `Tree{S}`; `NearestNeighbors` KDTree behind a
  `NearestNeighborIndex`; fold `SortedTupleSet` into the automaton work; remove `print`/`println` from RRT.
- Numerics: unify `bisection`(golden-section)/`newton_method`/`dbisection`; dedup the 4 LMI builders in
  `ellipsoidal_transitions.jl`.

### Phase 2 — System ⬜
Type the untyped wrappers (`SymbolicSystem` 20×`::Any`); `SystemApproximation` typed callables + threaded
`num_substeps`; static/dynamic controllers as **types**; consolidate the trajectory zoo to `Trajectory{T}`;
unify the 3 automata behind one CSR interface (relocated to Symbolic in Phase 5).

### Phase 3 — Mapping ⬜
Settle vocabulary (mapping = universe, state-set = subset; rename `get_state_domain`/`get_source_domain`);
make state-sets self-contained (stop threading `(set, mapping)`); type grids.

### Phase 4 — Problem ⬜
Immutable, fully-typed problems (`transition_cost`); uniform `discretize_problem` signature; plot via
`MS.stateset`; type `ap_semantics` (fix load order).

### Phase 5 — Symbolic ⬜
Slim the 12-param `SymbolicModelList`; unify `TimedHybridSymbolicModel`; own the automaton; keep the clean
execution-backend layer as the template.

### Phase 6 — Optim (centerpiece) ⬜
Adopt the Phase-0 base across all ~20 optimizers (delete the boilerplate + drifted wrappers); replace the
`isa`-chain with a `problem → solver` registry; decouple abstraction / synthesis / concretization into
testable stages; typed options; standardize field names; generalize `MOI_wrapper`; remove `∉` type piracy.

## Ecosystem adoption (interface-first)

| Concern | Package | Status |
| :-- | :-- | :-- |
| Sets + set algebra | LazySets | dep (partial) → finish migration |
| n-ball volume (Γ) | SpecialFunctions | transitive → promote |
| Sorted containers / heaps | DataStructures | dep (optim only) → use in Utils |
| Nearest neighbour (RRT) | NearestNeighbors | new |
| Small graphs (modes, PCLF) | Graphs | transitive → promote |
| Symbolic automaton (hot) | *custom CSR* | keep, behind interface |
| Systems | MathematicalSystems / HybridSystems | keep |
| Progress/logging | ProgressMeter | dep → replace `print` |

## Verification

- **Primary net:** the golden-output regression harness — after every phase, re-run `problems/*` +
  `docs/src/examples/*` and assert the synthesized abstraction/controller unchanged (within tol).
- **Per phase:** `Pkg.test()` (or targeted files / `--fast`); format (`JuliaFormatter`); Aqua; docs build;
  `@inferred`/JET on hot paths; BenchmarkTools on the abstraction build (no perf regression).

## Out of scope (v0.2)

Research-grade numerical kernels (SDP/LMI ellipsoidal transitions, PCLF) are cleaned and typed but not
re-derived. `BipedRobot/`, `control_server/`, `bench/` are updated only as needed to keep CI green.
