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

### Phase 1 — Utils (biggest ecosystem win) ✅ DONE (NearestNeighbors + `set_algebra.jl` rename deferred)
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
- ✅ **Collapse duplicate accessors — DONE.** One name per concept: `get_dim` (was `get_dim`+`get_dims`;
  neither exported, so the set (`UT`) and mapping (`MP`) methods coexist without unqualified-collision),
  `get_volume` (dropped the rectangle-only `volume` alias), `get_sublevel_set` (dropped `expand`; kept the
  `*` operator sugar), `affine_transformation` (dropped the unused ellipsoid-only `transform`),
  `is_intersecting` (was `is_intersection` on rectangles + `is_intersected` on ellipsoids). Ported all
  call sites in `src/`, `problems/*`, `test/*`, `docs`, `bench`. Root `utils/` case-study + `Depreciated`
  scripts intentionally left to lag (outside the v0.2 library surface). Caught a latent internal caller:
  `grid.jl::get_volume(::Grid)` used the dropped `volume` alias.
  Verified: rectangle/ellipsoid/lazy-set/grid units green, regression fingerprint unchanged (9705631),
  linearization 18/18, lazy-ellipsoid 13/13, uniform-ellipsoid 22/22, Aqua 11/11.
- ⏸ **NearestNeighbors — deferred.** The RRT distance is a custom non-`Metric` function over `Ellipsoid`
  states, so a KDTree would fall back to linear scan anyway. Revisit only if a Euclidean-metric NN index
  is needed elsewhere.
- ✅ **Numerics: scalar optimizers unified — DONE.** The three were misnamed/misplaced: `bisection` is
  actually golden-section, and `dbisection` (the only one with production callers) was hidden inside the
  *set* file `ellipsoid_inclusion.jl`; `bisection`/`newton_method` had no callers outside their own tests.
  Consolidated into `src/utils/optim/scalar_optimization.jl` as `golden_section_search` / `newton_method` /
  `derivative_bisection` — algorithms byte-identical, existing kwargs preserved (no caller breakage), each
  docstringed. Deleted `optim/{bisection,newton_method}.jl`; merged their tests into one
  `test/utils/optim/scalar_optimization.jl` (+ a first direct `derivative_bisection` test). Verified via
  scalar-opt unit + ellipsoid inclusion/intersection.
- ✅ **LMI builders: safe dedup done; deep dedup deferred.** Hoisted the 5 duplicated local
  `eye(n) = LA.diagm(ones(n))` closures to one module-level def in `ellipsoidal_transitions.jl` (identical
  matrices, zero call-site churn; `LA.I` can't substitute since `eye` also fills identity *blocks* inside
  matrix literals). Verified numerically identical: lazy-ellipsoid path cost `959.3348169927308` unchanged,
  both ellipsoid abstractions 13/13 + 22/22. **Deferred** the deeper block-PSD helper extraction: the 4
  builders are numerically-pinned research kernels (per-builder tolerances 1e-8/1e-10, J vs log-det
  objectives, P-as-variable vs parameter) feeding the exact regression fingerprint — factoring the PSD
  blocks risks silent numeric drift for low structural payoff, and the plan marks numeric kernels
  "cleaned/typed, not re-derived".

### Phase 2 — System 🔄 (edits done, gate running)
- ✅ **Typed wrappers** (`linearization.jl`): `SymbolicSystem` 20×`::Any` → 20 type params (constructed
  positionally in one place, `problems/non_linear.jl`, so zero call-site churn);
  `AffineApproximationDiscreteSystem{S, LT, F}` with a typed inner-constructor closure.
- ✅ **Kernel approximations**: all 7 silent empty-body interface stubs → error stubs; the 10 concrete
  approximation structs (`GrowthBound`, `Linearized`, `CenteredSimulation`, `RandomSimulation`,
  `OverApproximationMap` × discrete/continuous) lost their `Union{Nothing, …}` system fields (no caller
  ever passed `nothing`) and `::Function` fields → type params; the magic `num_substeps = 5` (4×
  hardcoded in `discretize`, plus `ngrowthbound`) is now one `ST.DEFAULT_NUM_SUBSTEPS` constant threaded
  as a kwarg through every `discretize`.
- ✅ **Controllers — two orthogonal axes, hierarchy + trait.** The audit's "AbstractStaticController/
  AbstractDynamicController as types" collides with reality: the discrete/continuous axis is load-bearing
  (~20 solver fields typed `ST.AbstractDiscreteController`/`…ContinuousController` — automaton-level vs
  concrete-level controller), and Julia has no multiple inheritance. Resolution: keep the
  discrete/continuous *hierarchy*, add the static/dynamic axis as the **`controller_kind` trait**
  (`StaticKind`/`DynamicKind` singletons, error stub — no silent fallthrough). `get_closed_loop_trajectory`
  now dispatches on the trait into two typed `_closed_loop` methods, replacing the runtime
  `initial_state === nothing` branch. Every concrete controller (8 in src + 2 BipedRobot in problems/)
  declares its kind; wrapper controllers delegate. Dead `state_domain`/`input_domain` controller methods
  removed (never called); protocol stubs now error with the offending type; argument names
  `controller_state`/`measurement`.
- ✅ **Trajectories**: deleted the never-used `HybridTrajectory`. Kept `DiscreteTrajectory`/
  `ContinuousTrajectory` (load-bearing for branch-and-bound/Bemporad-Morari) and `ClosedLoopTrajectory`
  (generators/certifiers) — the audit's "consolidate to `Trajectory{T}`" would rewrite those solvers for
  no user benefit; revisit in Phase 6 if their APIs merge. `wrap_coord`/`get_periodic_wrapper` moved to
  their single home `src/utils/periodic.jl` (+ scalar `wrap_value` shared with
  `Mapping.PeriodicGridMapping`, which keeps its precomputed-index fast path).
- ✅ **Automata**: dead `{N, M}` params dropped from `AbstractAutomatonList` and all 5 subtypes (three
  impls hardcoded `{3,3}`, ProductAutomaton/QuotientAutomaton used `{0,0}`); silent empty-body required
  methods → error stubs; non-idiomatic `NewXxxAutomatonList` free functions → `XxxAutomatonList(n, m)`
  constructors; `compute_post!` (used by docs Getting Started + unit tests but implemented only for
  Sorted) got a generic `append!(targetlist, post(...))` default so it works for every implementation;
  `FastIndexedAutomatonList` added to the automaton unit tests (was untested). CSR unification +
  relocation to Symbolic stays in Phase 5.

### Phase 3 — Mapping 🔄 (edits done, gated jointly with Phase 2)
- ✅ **Vocabulary settled — "domain" fully retired.** The Symbolic accessor surface returned
  `AbstractStateSet`s under the old Domain-module names, with a straight duplicate
  (`get_state_domain` ≡ `get_source_domain` ≡ `Xset`): renamed to `get_state_set` /
  `get_retained_set` / `get_input_set` across Symbolic + consumers, deleted the duplicate alias
  (which the rename would otherwise have turned into a self-recursive method that *overwrote* the
  error stub — same signature). Root `utils/` scripts lag, per policy.
- ✅ **Silent stubs → error stubs** in `SymbolicModel` (`get_state_mapping`/`get_input_mapping`/
  `get_state_set`/`get_retained_set`/`get_input_set`/`get_automaton`/`is_determinized`/`metadata` —
  `metadata` previously fell through to a cryptic `has_metadata(::Nothing)` MethodError) and in
  `Grid` (`get_origin`/`get_h`); `abstract_state_set.jl`'s bare `error("not implemented")` now
  names the offending type.
- ✅ **Typing**: `DeformedGrid{N, T, F, FI, AT}` (was `f::Function`, `fi::Function`,
  `A::Union{Nothing, Any}` — the docstring claimed `SMatrix`, the field said `Any`);
  `Set{Any}` → `Set{NTuple{N, Int}}` in the `SetMinus` pos-enumeration.
- ✅ **Hygiene**: `empty_states!(::ImplicitStateSet)` stores `empty_region(N)` directly (was the
  awkward `∅ \ ∅`).
- 📌 **`(set, mapping)` bundling → deliberately moved to Phase 5.** The bundle's only real
  consumer is `SymbolicModelList`, whose 12 type params get grouped in Phase 5 — bundling now
  would churn the same struct twice. The two-arg `MP` API stays until then.
- (Already fixed in Phase 0: broken `PeriodicGridMapping` constructors, orphaned
  `time_grid_mapping.jl`.)

### Phase 4 — Problem 🔄 (edits done, gated jointly with Phases 2–3)
- ✅ **All 6 `ProblemType`s immutable** (verified: no in-place field mutation anywhere in the repo).
- ✅ **`transition_cost::Any` → `TC` type param restored** on `OptimalControlProblem{S, XI, XT, XC, TC, T}`
  (the docstring already promised `TC`; no caller writes explicit type params, so zero churn).
- ✅ **Uniform `discretize_problem(problem, Δt::Float64; num_substeps)`** — fixed the
  `ReachAndStayProblem` keyword-`Δt` outlier; all defaults now come from `ST.DEFAULT_NUM_SUBSTEPS`.
- ✅ **`ap_semantics` typed**: moved `@enum INCL_MODE` + `invert_incl_mode` (was the camelCase
  `_invInclMode`) to `src/utils/incl_mode.jl` — Utils loads before both Problem and Mapping, which
  resolves the load-order cycle that had forced `Dict{Symbol, Any}`. `Mapping` aliases
  (`MP.INNER === UT.INNER`, …) so every existing `MP.INNER` spelling keeps working;
  `CoSafeLTLProblem.ap_semantics::Dict{Symbol, UT.INCL_MODE}` (a `Dict{Symbol, Any}` with
  INCL_MODE values still converts on construction). Stale `Dionysos.Domain.INNER` docstring fixed.
- ✅ **Plot recipes decoupled** from system internals: `problem.system.X` → `MS.stateset(problem.system)`
  (5 recipes; `Problem` now imports MathematicalSystems).
- ✅ `BisimulationQuotientProblem` docstring/type drift fixed (`{S,X,D,R,P,G}` → `{S, X, R}`).
- Kept: the `Infinity` sentinel; the `trajectory_success(::CoSafeLTLProblem)` placeholder
  (documented as such).

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
