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
- ✅ **`SpecialFunctions.gamma` for ellipsoid volume** — promoted `SpecialFunctions` to a direct dep;
  removed hand-rolled `gamma_half_integer_from_dim` (`src/utils/sets/ellipsoid.jl`). Verified.
- ✅ **Type-stable `Tree{S,A}`/`NodeT{S,A}`** (`src/utils/data_structures/tree.jl`) — bound the unbound
  `NodeT` UnionAll; `action::A` (defaults `Any`) instead of a bare `::Any` field; typed accumulators.
  Removed `print`/`println` from `src/utils/search/RRT.jl` (→ `@debug`); snake_case rename
  (`kNearestNeighbors→k_nearest_neighbors`, `get_nNodes→get_n_nodes`, `is_leave→is_leaf`, …). Verified
  (tree 68/68, lazy-ellipsoid 13/13).
- ⬜ **Sets → LazySets (staged; user chose full migration).** Target: delete
  `src/utils/sets/lazy_set_operations.jl` and represent set algebra with pure LazySets types.
  - **Correctness constraint (discovered):** the `A\B` discretization is inclusion-mode-aware —
    `states(A\B, incl) = states(A, incl) \ states(B, inv(incl))` (B uses the *inverted* mode). So a
    plain `Intersection(A, Complement(B))` discretized with one mode is WRONG. Canonical rep is
    `Intersection(A, Complement(B))` **and** the Mapping discretization must split it and invert B's
    mode. Also `set_in_period` (periodic wrapping) has no LazySets equivalent — it stays a Dionysos
    function that *produces* a `UnionSetArray`.
  - **Consumer contract to preserve:** construction in `problems/*` (`UT.LazySetUnion(v)`,
    `UT.LazySetMinus(A,B)` — path_planning, pwa_sys, articulated_vehicle, Pendulum/*); dispatch+fields
    in `src/mapping/grid_mapping/grid.jl` & `grid_mapping.jl` (`get_pos_from_set`/
    `get_states_from_set_strict` reach `.sets`, `.A`, `.B`); `src/mapping/abstract_state_set.jl`
    `ImplicitStateSet` (holds the minus, corner INNER/OUTER/CENTER via `.A`/`.B`, `set_in_period`).
    `HyperRectangle` doubles as an **integer index-box** (`get_pos_lims` returns `rectI` with Int
    bounds) — that use never enters a LazySets container, so it stays as-is for now (a later `IndexBox{N}`
    split is optional).
  - **Stages (each verified against fast set/mapping tests + regression net before the next):**
    1. ✅ `abstract type AbstractSetNode{N,T} <: LazySets.LazySet{T}` + `LazySets.dim`; old algebra kept
       working. Joining the `LazySet` hierarchy required constraining every leaf/container membership to
       `Base.in(x::AbstractVector, …)` (else it's ambiguous with LazySets' `in(::AbstractVector,
       ::LazySet)`) — done in `rectangle.jl`, `ellipsoid_inclusion.jl`, `lazy_set_operations.jl`.
       Verified: lazy_set_ops 45/45, rectangle/ellipsoid/abstract_state_set/grid green.
       Also fixed a **pre-existing latent bug**: `HyperRectangle(lb, ub)` with `length(lb)≠length(ub)`
       infinitely recursed → `StackOverflowError` (masked by a `@test_throws Exception`); now a clean
       `DimensionMismatch`.
    2.–4. ✅ **DONE.** `set_union(v) → UnionSetArray`, `set_minus(A,B) → Intersection(A, Complement(B))`;
       gutted `lazy_set_operations.jl` down to the LazySets-backed helpers + `AbstractSetNode` (kept the
       filename). Mapping discretization (`get_pos_from_set`/`get_states_from_set_strict`) dispatches on
       `SetUnion`/`SetMinus`/`EmptyRegion` and destructures the minus via `minus_included`/`minus_hole`
       (which invert B's inclusion mode — the correctness constraint above). `ImplicitStateSet` stores a
       `Region` and uses the total extractors. Ported `problems/*`, `MOI_wrapper`, the trajectory
       certifier, `format_input_set` (→ `IntersectionArray`), and all tests.
       **Gotchas hit & fixed:** (a) joining `LazySet` forced `Base.in(::AbstractVector, …)` on every
       leaf/container (ambiguity with LazySets' `in`); (b) LazySets' smart constructors simplify
       `Intersection(∅,·)→∅` and `Intersection(·,Universe)→·`, so `minus_included`/`minus_hole`/`add_set`/
       `remove_set` are **total** over any region, and `set_minus(A, ∅)=A` short-circuits (an empty
       `UnionSetArray` has no inferable `dim`, which otherwise trips `Intersection`'s dim assertion for
       obstacle-free domains e.g. `pwa_sys` simple mode); (c) `UnionSetArray` needs a `Vector{<:LazySet{T}}`
       not the unparametrised `LazySet`; (d) `UniformEllipsoidAbstraction` reached `X.A` → now
       `UT.minus_included(X)`; (e) the minus plot recipe needed a strictly-more-specific signature to beat
       LazySets' `Intersection` recipe (Aqua ambiguity).
       Verified: full `--fast` (34 files), regression net (fingerprint unchanged: 26285 / 9705631), Aqua
       11/11, lazy-ellipsoid 13/13, uniform-ellipsoid 22/22.
       *Optional follow-up:* rename `lazy_set_operations.jl` → `set_algebra.jl` (filename now stale).
- ⬜ Collapse duplicate accessors (`get_dim`/`get_dims`, `volume`/`get_volume`, `expand≡get_sublevel_set≡*`,
  `transform≡affine_transformation`, `is_intersection`/`is_intersected`).
- ⏸ **NearestNeighbors — deferred.** The RRT distance is a custom non-`Metric` function over `Ellipsoid`
  states, so a KDTree would fall back to linear scan anyway. Revisit only if a Euclidean-metric NN index
  is needed elsewhere.
- ⬜ Numerics: unify `bisection`(golden-section)/`newton_method`/`dbisection`; dedup the 4 LMI builders in
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
