# Plan — Multi-step walking as a hybrid system

Goal: turn the certified *single* footstep of `problems/BipedRobot/4D_model` into
**walking** — a certified infinite sequence of steps — by modelling the foot strike as
a discrete transition with a guard and a reset map.

---

## 1. The modelling insight: the state does not grow

Keep the state **role-relative** — `x = (hip_stance, knee_stance, hip_swing, knee_swing)`,
the stance foot pinned at the origin — and put the leg identity in the **mode**:

| | mode `L` (left foot planted) | mode `R` (right foot planted) |
| :--- | :--- | :--- |
| dynamics | `ẋ = u` | `ẋ = u` — **identical** |
| carved domain | `D` = "swing foot clear of ground and obstacles" | **the same `D`** |
| guard | swing foot in the ground band, ahead | the same |
| transition | `L → R`, reset `π` | `R → L`, reset `π` |

with the pair swap `π(x) = (x₃,x₄,x₁,x₂)`, obtained by rewriting the old swing chain from its
own foot once that foot becomes the new stance foot.

The point of this encoding: **the two modes are literally the same system**. Same dynamics,
same domain, same grid — because "the swing foot is clear" is the same predicate whichever
physical leg currently swings. So the abstraction is built **once** and used by both modes,
and the whole leg swap is carried by the reset map. Three consequences:

1. **No new dimension.** Neither the horizontal position `x` nor its derivative `ẋ` belongs
   in the state. `ẋ` has no meaning in a quasi-static, velocity-controlled model (there is
   no momentum for it to act on); `x` is a *derived* quantity (the real robot's boom
   position is reconstructed from the joint angles in `RSCore.x_to_mechanism_state`).
2. **The reset is exact on the abstraction.** Our grid has the same step on all four axes, so
   `π` maps cells bijectively onto cells: centres map to centres, `round` is exact, and the
   bisimulation survives the strike.
3. **The mode is the right home for the leg bit.** It is discrete, free for the synthesis,
   and it is what the driver and the real robot need to know.

The alternative — physical joint coordinates `(hipL, kneeL, hipR, kneeR)`, identity reset,
the permutation moving into the *domains* (`D` and `π(D)`) — is equally correct and closer to
the raw sensors, but the two modes then have different domains, so it needs a symmetry-aware
relabelling lift to avoid abstracting twice. The role-relative form gets the same reuse for
free.

---

## 2. The route: model it as a genuine hybrid system, and finish the machinery

The walking model *is* a hybrid system, and Dionysos already has the right concepts —
modes, guards, reset maps, per-mode grids. The machinery is simply **partial**, and the
gaps are small, well-identified, and useful beyond this problem. So the route is to build
the model the intended way and complete the library where it stops, rather than to route
around it with a bespoke lift.

What already works and needs nothing:

- self-transitions and multi-mode automata; per-mode dynamics, grids, time steps,
  `approx_mode`, `state_input_filter`;
- **per-mode carved domains** (`UT.SetMinus`) at the solver level — `build_mode_symbolic_models`
  passes each mode's own `X`, and holes are enumerated with the inverted inclusion mode,
  which is the sound choice;
- guards as boxes, half-spaces or arbitrary `LazySet`s, discretized `INNER` (exact for a
  `MP.CellUnion` guard);
- reachability, reach-avoid and safety on the flattened product automaton, which is an
  ordinary automaton the standard discrete solvers accept.

The missing pieces are listed in §5 and each is a self-contained improvement.

**One semantic point to settle rather than patch.** The abstraction turns a switch into an
*input the synthesis may decline*, while the front-end documentation describes it as taken
automatically inside the guard. For a foot strike, the controlled reading is the right one:
the impossibility of sinking through the ground is already enforced by the carved domain,
and *when* to put the foot down genuinely is a control decision in quasi-static walking. The
fix is to the documentation, not to the semantics.

---

## 3. The formal model

- **Continuous mode** (unchanged): `ẋ = u` on the carved domain `X ∖ infeasible_cells`,
  exact lattice (`τ·du = dx`), `CENTER_SIMULATION`, swept multi-cell steps.
- **Guard** `G` — the strike is available when the swing foot is on the ground and ahead of
  the stance foot:
  `G = cells_where(pos -> foot_y(centre) ≤ ground_band ∧ foot_x(centre) ≥ x_min)`.
  A `MP.CellUnion`, so its discretization is exact under any inclusion mode — no
  measure-zero manifold problem (an equality guard `foot_y = 0` would have no `INNER` cell
  at all; the band is the standard remedy and the ground carving already keeps that layer).
- **Reset** `π(θ) = (θ₃,θ₄,θ₁,θ₂)`, applied on grid *positions*, not coordinates:
  `π(p) = (p₃,p₄,p₁,p₂)`. Exact provided the four axes share the same step — asserted, not
  assumed.
- **Strike input**: one new abstract input symbol `σ_strike`, enabled **only** on `G`, with
  the deterministic transition `q → π(q)`.

**Why controlled switching is acceptable here.** The audit notes that Dionysos treats a
switch as an input the synthesis may decline. Physically the foot *cannot* sink through the
ground — but that constraint is already enforced by the ground carving. What remains is
"when do I put the foot down", which in quasi-static walking genuinely *is* a control
decision. So the controlled-switch semantics is the right one for this model, not a gap to
patch.

---

## 4. The specification: a gait is a *recurrence*, not an invariant

Chaining `N` single-step problems certifies `N` steps. The right infinite statement is that
the strike can always be taken again.

**Safety is the wrong specification, and cheaply so**: `u = 0` is an admissible input, so
staying inside the domain forever is satisfied by standing still. The largest
controlled-invariant set is essentially the whole domain, and certifies nothing about
walking. Reach-and-stay fails for the same reason.

Recurrence — "strike infinitely often" — is a Büchi condition, and the discrete solvers cover
co-safe LTL (finite words), not full LTL. But a **finite** computation certifies it. Let `G`
be the guard and `Reach(G)` the set of states from which `G` is reachable. Then

> if `π(G) ⊆ Reach(G)`, the robot can strike forever

by induction: from a post-strike state reach the guard again, strike, repeat. That is one
reachability synthesis on the hybrid product, with the solvers already wired in and no Büchi
machinery. The controller is the reachability controller re-targeted after each strike —
exactly the shape of the 6D model's two-step walking controller.

This also exposes, rather than hides, the classical sequential-composition condition
(Burridge–Rizzi–Koditschek): the post-strike image must land in the basin of the next swing.
When the inclusion fails, the check *names* the post-strike states that are not recurrent,
and the remedy is a wider guard or a finer grid — not a change of model.

**Verified** on the coarse model (`dx = 0.1`, `|θ| ≤ 0.6`, no obstacle): the synthesis
succeeds and every post-strike state is recurrent, in both stance modes.

---

## 5. The missing pieces, and the phases that add them

Each phase is independently useful and independently testable. Sizes are rough.

### Phase 0 — two soundness/robustness fixes (small, no new feature)

| Piece | Where | Why |
| :--- | :--- | :--- |
| flat automaton sized by the **alphabet**, not by the count of *used* inputs | `symbolic/lifts/hybrid_transition_assembly.jl` (`ninputs = length(inputs_set)`) | a `state_input_filter` leaves gaps in the used-id set, so `max(used_id) > nsymbols` and the discrete solvers index a dense `nstates × nsymbols` table out of bounds. We *will* use such a filter. One line, plus a regression test. |
| assert the reset is **lattice-exact**, else refuse | same file, next to the reset application | the reset is applied to the *cell centre* and rounded to the nearest target cell. That is exact iff the reset is a lattice automorphism (identity and pair-swap are, when the permuted axes share a step). Off-lattice by half a cell you silently get a plausible-looking, wrong abstraction. Mirror `MP.is_lattice_exact`. |

Optional follow-up, larger: a genuinely **over-approximated reset** (image box →
`get_states_from_set(..., MP.OUTER)`, several transitions per source cell) so that nonlinear
or off-lattice resets become sound instead of forbidden. The plumbing already accepts several
transitions per source.

### Phase 1 — modes that share one abstraction

Today `build_mode_symbolic_models` is an unconditional `map`: every mode gets a full
`MOI.optimize!`, even when two modes are *the same object* — the existing test passes
`[mode, mode]` and abstracts it twice. With `N` modes over `M < N` distinct dynamics you pay
`N` abstractions. For the walking model, whose two modes are the same system, that is exactly
a factor two on the dominant cost.

The feature: **let the user declare which modes share an abstraction**, e.g. a per-mode
`shared_abstraction` entry (`[nothing, 1]` = "mode 2 reuses mode 1's model"), resolved in
topological order before the build. Explicit beats implicit here — the user knows two modes
are the same system, and a declaration is checkable and self-documenting.

Worth adding underneath as a safety net: memoize on `(objectid(physical system), kwargs)`, so
that passing literally the same system with the same configuration is never abstracted twice
even without a declaration. That alone fixes the `[mode, mode]` case.

Validation when a declaration is made: the two modes' systems must agree on dynamics, state
set and grid configuration — otherwise refuse, rather than silently abstract the wrong thing.

Tests: a 2-mode system whose modes share one abstraction builds the underlying abstraction
once, and its flattened automaton is identical to the one obtained without sharing.

**Known residual cost, to measure rather than assume.** Sharing removes the duplicated
*build*, not the duplicated *product*: `add_intra_mode_transitions!` still copies each mode's
transitions into the flat automaton, so the flat automaton stays ~2× at ~40 M transitions per
mode. If that proves to be the binding constraint, the options are a lazy/shared product
representation, or the single-mode variant (one mode, one self-transition, leg bit tracked by
the driver) which avoids the product entirely — at the cost of losing the declarative mode
structure and hitting the parallel-transition label bug if two strike variants are ever
needed.

### Phase 2 — the walking model (direct MOI)

Two modes `L`/`R`, identical dynamics, domains `D` and `π(D)`, guards as `MP.CellUnion`s
(swing foot in the ground band and ahead), **identity** resets, mode `R` derived from mode
`L` by `PermutationLift`. Specification: `SafetyProblem` for the gait invariant of §4 — it is
already supported end to end, so this phase needs no further library work.

Test: the invariant set is non-empty, contains the nominal posture, and a closed loop of ≥ 5
strikes stays in it.

### Phase 3 — make it expressible in JuMP (finish the wrapper)

The front-end is partial on the hybrid path; three gaps, in order of importance:

| Gap | Where | Fix |
| :--- | :--- | :--- |
| `x ∉ O` **silently dropped** for hybrid models — a mode's domain is built from its box bounds only | `wrapper/lower_hybrid.jl` `_build_mode_system`; `obstacle_sets` is only called on the monolithic path | carve per mode, and allow the `∉` constraint to be mode-scoped so each mode gets its own obstacles |
| specification vocabulary limited to `Final` / `Always` | `build_hybrid_problem` | accept `EventuallyAlways` (reach-and-stay) once Phase 4 wires the solver |
| initial state must be a **single point** | `_hybrid_initial_state` | accept a region, as the monolithic path does — needed to state a basin or an invariant set as the initial condition |

Then the model of Phase 2 is written declaratively: `@mode`, `add_transition!` with a
`Guard(...)`, identity reset, `set_attribute(mode, "state_grid", …)`.

### Phase 4 — reach-and-stay on hybrid models

`control_solver_for` accepts only `OptimalControlProblem` and `SafetyProblem`. The flattened
automaton is an ordinary automaton, so `OPDS.OptimizerReachAndStayProblem` runs on it
already; what is missing is two lines of wiring **and** a dynamic controller concretizer —
today `build_concrete_controller` always builds the static `HybridQuantizedStaticController`,
which cannot represent a controller whose output depends on a specification-automaton state.

### Phase 5 — driver, and obstacles at absolute positions

`examples/BipedRobot/biped_4d_walking.jl`: closed loop over several strikes, accumulating the
world position driver-side (not a state) so the robot visibly advances; the dashboard's
`state_background!` hook already draws the carved slice.

Then the genuinely position-dependent case: with a fixed obstacle the carved domain depends
on where the stance foot is relative to it, so the world position enters as a **discrete**
index — a short chain of modes (far, approach ×k, straddle, past), each with its own carved
domain, chained by the same identity reset. Only the modes whose swing workspace actually
meets the obstacle need their own abstraction; the rest are permutation-derived.

Cheap alternative to evaluate first: a **step-primitive library** — the abstraction does not
depend on the target, so one abstraction plus one cheap synthesis per relative foothold Δ
gives a certified catalogue of steps for a footstep planner to chain. A Δ with a vertical
component gives stairs and uneven terrain for free, since `target_set` already accepts an
arbitrary Cartesian foothold.

---

## 6. Known traps (from the audits)

| Trap | Where | Status |
| :--- | :--- | :--- |
| Hybrid front-end silently drops `x ∉ O` | `wrapper/lower_hybrid.jl` | fixed in Phase 3; until then use direct MOI, which carves correctly |
| Reset is centre-point, never over-approximated | `hybrid_transition_assembly.jl` | identity/permutation are exact on our grid — asserted in Phase 0 |
| Flat automaton sized by *used* inputs | `hybrid_transition_assembly.jl` | fixed in Phase 0 (we trigger it via `state_input_filter`) |
| Parallel transitions share a label `"SWITCH s -> t"` | `hybrid_system_abstraction.jl` `_find_switch_transition` | avoided: `L→R` and `R→L` are distinct pairs. Would bite the one-mode variant |
| Reset image outside the domain ⇒ edge dropped silently | `hybrid_transition_assembly.jl` | count and report; a fully dropped guard must be an error, not a warning |
| Hybrid initial set must be a single point | `lower_hybrid.jl`, `optimal_control_problem.jl` | Phase 3 |
| Flat state ids depend on `Set` iteration order | `flat_index.jl` | never key a test or artifact on flat ids |
| Switch documented autonomous, implemented controlled | `wrapper/modes.jl` vs `global_input_map.jl` | documentation fix — the controlled reading is the correct one here (§2) |

---

## 7. What this does *not* do

No dynamic walking: no momentum, no impact law, no balance criterion (ZMP, capture point).
The quasi-static assumption is what makes the kinematic model legitimate, and it is the
assumption to revisit — with the torque models — before claiming anything about a fast gait.
The certified statement here is: *the joint trajectory reaches successive footholds while
keeping the swing foot clear of ground and obstacles, at sampling instants and in between,
with a trackable velocity profile.*
