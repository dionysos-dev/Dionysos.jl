# Plan — refactoring `PCLFBisimulationQuotient` into a modular solver

> Scope: `src/optim/hybrid_systems/PCLFBisimulationQuotient/` (1 905 lines, 83 definitions) plus the
> two `Utils` files it leans on. **Purpose: structure and clarity, not speed** — the performance
> work is finished and §6 explains why there is little left to win.
>
> **Status — P1 to P4 executed (2026-08-26), P5 deliberately not.** `atol` stays at `1e-4`.
> Three of the findings below were wrong and §8 records what executing them actually showed;
> they are left in place because the corrections are the useful part.
>
> This is a separate document because `plan.md` is the live audit of `HybridSystemAbstraction`
> (issue #501) and still has P5–P8 outstanding. Merge the two, or supersede that file, only
> deliberately.

## 1. Where the module stands

Nine commits on `jc_pclf_abstraction` fixed two real defects, one calibration and one file. What
remains is structural.

| file | lines | defs | state |
| :--- | ---: | ---: | :--- |
| `pclf_bisimulation_quotient.jl` | 631 | 15 | **cleaned** — documented, sweep flattened, names fixed |
| `bisimulation_quotient.jl` | 521 | 33 | three concerns in one file; six-fold accessor duplication |
| `cosafe_ltl_problem.jl` | 484 | 21 | four parallel lookups; `_lifted` naming with no counterpart |
| `sublevel_support.jl` | 269 | 14 | dispatch avoided in favour of name-mangling |

Already fixed and not revisited here: the quotient-volume algorithm (≈1000× and de-biased), the
`success` flag, the `atol` calibration, the module entry file name, `ST.mode_matrices`, and the
construction algorithm's documentation, nesting and naming.

## 2. What is actually wrong

**F1 — `bisimulation_quotient.jl` holds three unrelated things.** The data structure
(`PCAbstractState`, `PCBisimulationQuotient`, `add_state!`, `remove_state!`, `add_transition!`),
the volume geometry, the plot recipes, and **22 of its 33 definitions are statistics**
(`states_by_*`, `outgoing_degree_*`, `deadend_states`, `cell_complexities`, `bisimulation_stats`,
`print_bisimulation_stats`, …). Two thirds of the file is reporting sitting on top of the type
everything else in the module depends on.

**F2 — six accessors that are two.** `states_in_slice` / `states_in_obs` / `states_in_node` and
`state_ids_in_slice` / `state_ids_in_obs` / `state_ids_in_node` each filter on one field. The
`state_ids_*` trio additionally carries a `haskey` guard the `states_*` trio omits, so the two
families disagree about whether stale ids can reach them — a real inconsistency, not just repetition.

**F3 — most of that surface is unused.** Measured across `src/`, `test/`, `research/` and
`examples/`: `states_in_slice`, `states_in_obs`, `states_in_node`, `state_ids_in_slice`,
`state_ids_in_obs` and `num_faces` have **zero** callers outside their own file;
`state_ids_in_node` has 3 and `cell_complexities` has 2. Most of it is API that was never wanted.

**F4 — `sublevel_support.jl` encodes types in names.** `_support_abs_row_on_hyperrectangle` and
`_support_abs_row_on_hpolytope` differ only in an argument type that is *already in the signature*.
`CLAUDE.md` names multiple dispatch as the house idiom; here the caller picks the method by hand, so
supporting a third set type means editing every call site instead of none.
`radius_hyperrectangle(X) = LazySets.radius_hyperrectangle(X)` is a wrapper that adds nothing but a
shadowing name.

**F5 — `cosafe_ltl_problem.jl` naming.** `solve_concrete_problem` / `solve_concrete_problem_lifted`
and `initial_controller_memory` / `initial_lifted_controller_memory` are *not* duplicated logic —
the short name unpacks the optimizer and delegates. But `_lifted` implies a non-lifted counterpart
that does not exist, and four near-parallel lookups (`_find_qid_in_node`, `_find_qid_same_node`,
`_find_qid_global`, `_find_successor_qid`) want reading before anyone trusts them.

**F6 — an unresolved author note ships in the code.** `# maybe something to change here??` at
`pclf_bisimulation_quotient.jl:337`, on the `set_difference_decompose` call inside
`refine_state_by_observation!`. Left deliberately: it flags the `atol` question §5 records.

## 3. Phases

Ordered so each is independently revertible and the risky one comes last.

**P1 — split `bisimulation_quotient.jl`** (addresses F1). `git mv` the statistics block into
`quotient_stats.jl` and the plot recipes into `quotient_recipes.jl`, leaving the data structure and
volume geometry. Pure `MOV`; no line changes. Verify by state count and the existing suites.

**P2 — collapse the accessors** (F2, F3). Two functions with a field selector replace six. Delete
the six with no callers rather than port them — they are unused surface, and `git` keeps them.
Settle the `haskey` disagreement explicitly: decide whether stale ids are reachable and apply one
answer. This is the only phase that could change behaviour, hence the explicit decision.

**P3 — dispatch in `sublevel_support.jl`** (F4). One `_support_abs_row` with two methods; drop the
`radius_hyperrectangle` wrapper. Mechanical.

**P4 — read `cosafe_ltl_problem.jl` and report** (F5). Deliberately *not* a change phase. The four
lookups need understanding before touching; rename `_lifted` only once it is clear what the lift is.

**P5 — degeneracy filter in `set_difference_decompose`** (see §5). Optional, and the only phase with
a real chance of surprises.

## 4. What not to do

Four optimizations were implemented, measured and rejected this session. Recording them so they are
not attempted again:

- **Spatial index on the volume computation** — 1.1× over a tuned linear scan. Not worth it.
- **Bounding-box screen inside `poly_intersection_parts`** — *counterproductive*.
  `box_approximation` on an `HPolytope` costs LPs, and there each box faces only a few counterparts,
  so it never amortizes.
- **Hoisting `source_ids` out of the target loop** — not semantics-preserving. Refinement removes and
  adds states mid-loop, and the rescan is how freshly-split pieces become visible to later targets.
- **State-level bounding-box screen** — implemented, provably sound, **95.62 % skip rate, zero
  speedup**. This is the most useful negative result: the ~695 000 wasted refinement calls are
  *cheap* early-exits. Cost lives in the ~12 700 calls that actually fire, at roughly 4 ms each.

## 5. The `atol` question, unresolved

`atol` is the inset applied when one polytope is cut from another. It is bounded on both sides and
neither bound is comfortable:

| `atol` | domain left uncovered | outcome |
| :--- | ---: | :--- |
| `1e-3` (former default) | 3.71 % | silently drops a twenty-fifth of the space |
| `1e-4` (current) | 0.386 % | works end to end |
| `1e-6` | 0.004 % | **389 flat pieces survive; 4 have no support vector — plotting dies** |
| `0` | — | **diverges** (killed at 8 h / 4.7 GB) |

The loss compounds: measured on one run, the slice-level gap was amplified ≈15× by the time
refinement finished, because every split cuts again.

P5 would remove the ceiling: `set_difference_decompose` keeps any piece where `!isempty(piece)`, but
`isempty` is a *feasibility* test that accepts flat, measure-zero polytopes. Dropping pieces with no
interior costs nothing real (they contribute zero volume) and would make `1e-6` viable. The cost is
an interior test per surviving piece, and the risk is that "no interior" is a numerical judgement.

## 6. Why performance is not in scope

Post-cleanup profile of a 7-level, 2-mode construction: `refine_one_state!` is 74 % of the time,
split 41 % intersection / 24 % set-difference. Both are genuine geometry — the intersection half
already routes through LazySets' LP-free 2-D `HPolygon` path. With the state-level screen proving
attempted work is nearly free, the remaining levers are making each *firing* refinement cheaper or
producing fewer of them, and the latter means a coarser quotient. Parallelism is the only avenue
that could give a multiple rather than a percentage; `src/symbolic/` has threaded and distributed
backends, this module has none, and the sweep mutates shared state so it is not a small change.

## 7. Verification

The unit suites are necessary but weak. The sharp check is that the example's construction still
yields **10 611 states** — a restructure that subtly reordered refinement shows up there and nowhere
else. Run per phase:

```
julia --project=test test/optim/PCLFBisimulationQuotient/bisimulation_quotient.jl
julia --project=test test/optim/PCLFBisimulationQuotient/recipes.jl
julia --project=test research/BisimulationQuotient/paper_example_3_1.jl  # states, volume, plots
```

The example matters because it is the only thing exercising the plotting path, which is exactly
where `atol = 1e-6` failed while construction and unit tests passed.

Current reference values: 10 611 states, 0 deadends, volume **186.72811640312693**,
`Initial set controllable: true`, construction ≈88 s.


## 8. What executing it showed

**F2 was wrong.** The `haskey` asymmetry between `states_in_*` and `state_ids_in_*` is correct by
construction, not an inconsistency: every caller passes `state_ids = controllable_set`, a
caller-supplied subset that genuinely can name removed states, whereas `values(T.states)` cannot.
P2 became a plain deletion of five accessors with no callers anywhere.

**F4 was wrong about the cause.** `_support_abs_row_on_hyperrectangle` was not a dispatch
mistake but dead code, reachable only from a commented-out line, with the local
`radius_hyperrectangle` wrapper existing solely to serve it. Deleting would also have been wrong:
that comment sits in a function whose `X` is typed `HPolytope`, so it is a remnant of when `X`
could be a box. Two methods of one `_support_abs_row` keep the closed-form path and let dispatch
choose — which is what the comment was reaching for.

**F5 was wrong.** The four lookups are not near-duplicates but a cost-ordered fallback cascade:
successors of the current state, then the rest of its node, then the whole quotient. And the
reason it exists is the `atol` gaps of §5 — a trajectory can land in a gap belonging to no state,
or drift across a boundary into a state that is not a predicted successor. The wide searches are
what keep simulation running, so the gaps and the cascade are one story.

**Two things the plan missed entirely.**

- `print_bisimulation_stats` (5 callers) recomputed every figure that the dead
  `bisimulation_stats` already assembled. Deleting the dead one would have entrenched the
  duplication; the printer now renders from it. Output verified byte-identical.
- `simulate_closed_loop` re-implemented mode-matrix extraction inline, on every step, duplicating
  `ST.mode_matrices`. It now extracts once.

**A method note.** The first audit counted callers outside each file and reported 14 unused
functions; counting internal callers too left 4, and one of those was the better half of a
duplicated pair. *Unused* and *dead* are not the same question, and only the second justifies
deletion.

**Still open.** `cosafe_ltl_problem.jl` and `sublevel_support.jl` are read only in the places
touched here — roughly half of each remains unexamined. `part_ids` is a `Vector` used as a set
(`remove_state!` filters it linearly; profiled under 1%, so wrong container rather than hot).
`solve_concrete_problem_lifted` is explicable now — the lift is the quotient's sparse ids to a
dense `1:n` automaton — but still misnamed.