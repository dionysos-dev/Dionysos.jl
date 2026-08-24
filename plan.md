# Plan — audit and generalization of `HybridSystemAbstraction`

> Addresses issue **#501**. Supersedes the Adaptive Cruise Control plan that lived here; that
> work is merged and the old document is recoverable at `git show 9e729350f:plan.md`.

> **Status — P0 through P4 are built** (2026-08-24), which is the whole programme §5 recommends.
> F1, F2 (documentation half), F3, F5, F9a, F9b, F9c, F9e, F9f and F9g are closed; hybrid
> reach-and-stay works end to end; the initial condition is a set; a switch that discretizes to
> nothing now fails instead of silently disconnecting the modes; and the modes can be abstracted
> in parallel. **P5–P8 remain, all gated on a driving example** (§5, decision 3).
>
> Building it turned up one thing the audit had not: the empty-abstract-set check added for F9g
> fired immediately on the repository's own test suite. `test/wrapper/hybrid.jl`'s half-space-guard
> test asked for a `[1.0, 1.5]²` target on a `0.5` grid, which under `INNER` inclusion needs a cell
> centred at `1.25` — not a grid point. That test had been solving against an **empty target set**
> and passing, because it only asserted that the abstraction was built. Widened, and it is now a
> real end-to-end check.

**Scope.** The whole hybrid path, end to end: the JuMP front-end that lowers a moded model
(`src/wrapper/{modes,parse_scoped,clock,lower_hybrid,solver_selection}.jl`), the abstraction it
builds (`src/symbolic/lifts/`), the solver family
(`src/optim/hybrid_systems/HybridSystemAbstraction/`), and the closed-loop simulation.

**How the findings were established.** Everything below is read from the code, and the three
findings marked **[verified]** were additionally reproduced by running the solver. Findings
marked **[read]** follow from an unambiguous code path but were not executed. Nothing here was
measured: the performance claims (F7, F8) are arithmetic on the data structures, not benchmarks.

**Where the defects are.** Five of the eight findings are in the *lowering* — `src/wrapper/`, not
the solver. The solver core comes out of this well: the lift architecture composes cleanly, the
mode-sharing mechanism is careful, and the reset-quantization warning shows the authors had
already found the sharpest edge in it. Effort belongs in the front-end.

---

## 1. What the solver does today

```
JuMP model with @mode/add_transition!
  → ModelIR                                       (wrapper/model_ir.jl)
  → HybridSystems.HybridSystem + HybridSpec       (wrapper/lower_hybrid.jl)
  → one UniformGridAbstraction per mode           (hybrid_symbolic_builder.jl)
  → ClockLift per mode, then ModeLift             (symbolic/lifts/)
  → HybridSymbolicModel  (one flat automaton, one global input alphabet)
  → a discrete solver on that automaton           (optim/discrete_systems/)
  → HybridQuantizedStaticController
```

The two structural ideas are worth naming, because the whole design rests on them and both are
sound:

- **A mode is just a system.** Each mode is abstracted by the ordinary
  `UniformGridAbstraction.Optimizer`, so every feature of the continuous solver — growth bounds,
  approximation modes, parallel build backends — is inherited for free, per mode.
- **The composition is a lift.** `ClockLift` (add a time axis) and `ModeLift` (glue modes through
  guarded switches) are `AbstractLift`s producing another `AbstractSymbolicModel`, so the result
  is fed to the *same* discrete solvers as a flat abstraction. Nothing about reachability or
  safety synthesis is re-implemented for the hybrid case.

That is a good architecture. Most of what follows is about the seams.

---

## 2. Capability inventory

| | Supported today | Where it is decided |
| :--- | :--- | :--- |
| Reach (`Final`) | ✅ | `solver_selection.jl:45` |
| Reach-avoid (`Final` + `Always`) | ✅ | `lower_hybrid.jl:386` |
| Safety (`Always`) | ✅ | `lower_hybrid.jl:398` |
| Reach-and-stay (`EventuallyAlways`) | ✅ *(was: silently dropped)* | F1, P2 |
| Co-safe LTL over modes | ❌ rejected | `solver_selection.jl:45` |
| Abstraction only (no spec) | ❌ from JuMP, ✅ from MOI | `lower_hybrid.jl:401` |
| Controller-chosen switch | ✅ (but documented as unsupported) | F2 |
| Environment-forced switch | ❌ **not expressible** | F2 |
| Non-box guards | ✅ time-free only | `lower_hybrid.jl:171` |
| Reset maps | ⚠️ sound only if lattice-exact | F4 |
| Initial **set** | ✅ *(was: collapsed to a point)* | F5, P1 |
| Mode with no `Always` | ✅ unconstrained *(was: forbidden)* | F3, P1 |
| Parallel mode abstraction | ✅ opt-in `parallel_modes` | F9b, P4 |
| Clocks | ⚠️ exactly one, all modes or none | F6 |
| Per-mode grids / time steps | ✅ | `optimizer.jl:223` |
| Shared abstraction between modes | ✅ | `hybrid_symbolic_builder.jl:64` |

---

## 3. Findings

Ordered by severity. Severity = (how wrong the answer can be) × (how likely the user hits it
without noticing).

### F1 — `EventuallyAlways` on a mode is parsed, then silently discarded **[verified]**

`parse_scoped.jl:212-239` accepts any `SpecSet` on a mode, and `SpecSet` includes
`EventuallyAlways` (`specification.jl:87`). The entry lands in `mode.specs` with kind
`EVENTUALLY_ALWAYS`. But `build_hybrid_problem` (`lower_hybrid.jl:378-405`) only ever reads
`FINAL` and `ALWAYS`. The entry is dropped without a word.

Reproduced — the same model lowers identically with and without the constraint:

```
without EventuallyAlways : SafetyProblem
with    EventuallyAlways : SafetyProblem
  -> identical problem type? true
```

So a user who writes `Always(S) ∧ EventuallyAlways(T)` on a hybrid model gets a controller for
`Always(S)` alone, reported `OPTIMAL`, with no indication that half the specification was
ignored. Written alone it is worse: the model falls through to
`"A hybrid model needs a specification on at least one mode"`, which names the one thing the user
did do.

This is the issue's first bullet, and it is not merely "not implemented" — it is *accepted and
ignored*, which is the failure mode a solver library must never have.

**Fix.** Two parts, and the first must land regardless of whether the second ever does: reject
unhandled spec kinds explicitly in `build_hybrid_problem` (name the kind and the mode), then
implement reach-and-stay properly (P2 below — it is cheap).

### F2 — The switch is controller-chosen, but documented as environment-driven **[verified]**

`modes.jl:284-303` rejects `ControlledSwitching` and accepts only `AutonomousSwitching`, with
this rationale in the error and again in `README.md:320`:

> "The hybrid abstraction takes a switch whenever the state is in the guard, so a controlled
> switch would silently behave as an autonomous one."

The implementation does the opposite. `GlobalInputMap` (`global_input_map.jl:50-64`) allocates
**one input symbol per hybrid transition**, and `add_inter_mode_transitions!`
(`hybrid_transition_assembly.jl:53-108`) labels every switch with that symbol. Meanwhile
`add_intra_mode_transitions!` (`:23-42`) copies each mode's own transitions in **unconditionally
— including for states inside a guard**. So from a guard state the automaton offers both "switch"
and "stay under some continuous input", and the discrete solver picks one. The controller decides.

Reproduced on a two-mode model whose guard covers the source mode's entire state set:

```
total inputs = 7  (continuous 1:6, switching 7:7)
mode-1 states total                        : 7
  with an intra-mode successor (may stay)  : 7
  with a switch successor                  : 7
  offering BOTH => controller may decline  : 7
```

Every state that fully satisfies the guard can still decline the switch.

**This is a documentation inversion plus a missing feature — not a defect in what is
implemented.** The may-semantics the code implements is correct, and it is what every model in
this repository actually wants. The clearest evidence is the 4-D biped, whose own docstring
([`walking.jl:72-74`](problems/BipedRobot/4D_model/walking.jl)) says:

> "The switch is offered to the synthesis rather than forced: the impossibility of sinking
> through the ground is already carried by the carved domain, and in quasi-static walking *when*
> to put the foot down is genuinely a control decision."

— while passing `AutonomousSwitching()` to `HybridSystems` on line 85 as a placeholder. So the
repository's most sophisticated hybrid model documents the opposite of what `modes.jl` claims,
and relies on the behaviour `modes.jl` says is unsupported.

The cost is borne by a *future* user: someone modelling a fault, a bouncing ball, or a threshold
that physically forces the mode change reads `AutonomousSwitching` in the API, believes the switch
is forced, and gets a controller that assumes it can refuse. Nothing flags it. But no such model
exists here today, which is what demotes the fix from urgent to eventual.

**Fix.** P0 (documentation, immediate) then P7 (implement forced switching) — *if* a model that
needs it ever appears. See §5.

### F2b — There is no notion of a mode invariant **[read]**

A hybrid automaton in the standard formulation pairs each transition's **guard** with each mode's
**invariant**: you may occupy mode `k` only while `Iₖ` holds, and it is the guard/invariant pair —
not the guard alone — that makes a switch forced. Dionysos has no invariant concept. A mode has a
state set (its box bounds, `_mode_box` at `lower_hybrid.jl:17`), and a trajectory that leaves it
simply produces an abstract state with no outgoing transitions: a **deadlock**, not a forced
switch.

This is more fundamental than the may/must distinction of F2, and it constrains any future
attempt at forced switching: "the switch must be taken before the invariant is violated" is not
expressible, so an urgent-transition semantics would have to invent the invariant first. Noting
it here so P7 is not designed around the guard alone.

### F3 — A mode with no `Always` becomes entirely unsafe **[verified]**

`_hybrid_spec` (`lower_hybrid.jl:327-343`) skips any mode with no set of the requested kind, and
`states_satisfying(::HybridSymbolicModel, ::HybridSpec)` (`states_satisfying.jl:26-36`) only
contributes states for modes present in `per_mode`. `satisfies` agrees
(`specifications.jl:99-103`: absent mode → `false`).

Reproduced — `Always` written on mode `a` only:

```
safe spec covers modes   = [1]
satisfies(x=2.0, mode 1) = true
satisfies(x=2.0, mode 2) = false
```

For a **target** set that reading is right: a mode you said nothing about is not a goal. For a
**safe** set it is backwards: saying nothing about mode 2 becomes "mode 2 is forbidden entirely",
so the synthesis quietly refuses to ever enter it. On a model with more than two or three modes
this is very easy to hit and presents as unexplained infeasibility.

**Fix.** Make the default per kind, not per spec: absent mode ⇒ *empty* for `FINAL`, *the mode's
whole state set* for `ALWAYS`. One branch in `_hybrid_spec`, plus the matching default in
`satisfies`.

### F4 — Reset maps are sound only when they are lattice-exact **[read]**

`add_inter_mode_transitions!` applies the reset to the source cell's **centre**
(`hybrid_transition_assembly.jl:91-93`) and quantizes the image to a **single** target cell. The
true image of a *cell* under a non-affine, or affine-but-unaligned, reset is a set that may span
many target cells; recording one of them is an under-approximation of the successor set, which is
exactly what soundness forbids.

The code knows this. `_reset_quantization_offset` / `_warn_reset_not_lattice_exact`
(`:128-155`) measure the snap distance and warn, with an accurate diagnosis in the message
("the abstraction may be unsound"). But a `@warn` with `maxlog = 1` in a long solver log is not a
guarantee, and there is no way to ask for the sound behaviour.

**Caveat on severity — no model in this repository is exposed.** The thermostat's reset is the
identity; the biped's is `leg_swap(θ) = SVector(θ[3], θ[4], θ[1], θ[2])`, a coordinate permutation
on a uniform `dx = 0.1` grid — exactly the lattice-exact case the warning exempts. So this is a
landmine for a future user, not a present defect, and it has **not** been demonstrated: no
counterexample model was constructed. Building one is the first task of P5, before any fix.

**Fix.** Over-approximate the reset image of the whole cell and emit a transition to every
intersecting target cell (P5). The machinery already exists on the continuous side —
`System.approximation/` computes reachable-set over-approximations from a growth bound, and a
reset map is the degenerate zero-time case.

### F5 — The initial condition collapses to a single point **[read]**

`_hybrid_initial_state` (`lower_hybrid.jl:357-370`) takes `LazySets.center` of the declared start
box, and both hybrid sub-solvers then take a single abstract state
(`optimal_control_problem.jl:136-139`, `safety_problem.jl:93-96`: `abstract_initial_set = [q0]`).

So on a hybrid model `Start(S)` means "start at the centre of `S`", `TerminationStatus` answers
for that one point, and a controller reported `OPTIMAL` may fail from a neighbouring state inside
the very set the user declared. The continuous solver has no such restriction. The comment at
`:353-356` states the design intent ("a hybrid problem starts from a point rather than a region")
but nothing enforces or surfaces it at the API level.

**Fix.** Accept a set (P1): `states_satisfying` already returns a state *list* for a `HybridSpec`,
and the discrete solvers already take an initial *set*. Not merely deleting the collapse, though —
the augmented initial state is also what `_simulate_hybrid` starts from (`results.jl:90`) and what
`trajectory_success` checks against, so the point and the set have to be separated rather than
swapped.

### F6 — One clock, all modes or none, and mixed models are structurally impossible **[read]**

Three separate restrictions compound (the issue's bullets 2 and 4):

- `detect_clock!` (`clock.jl:66-108`) errors on more than one clock candidate.
- `_build_mode_system` (`lower_hybrid.jl:74-76`) consults the *global* `clock_index`, so a clock
  is added to **every** mode, and `clock_system` (`clock.jl:126-136`) then errors on any mode that
  did not declare a rate.
- The arity mismatch is real, not incidental: `_reset_function` (`lower_hybrid.jl:214-239`) is
  built once per transition against the global clock, and `get_next_aug_state`
  (`hybrid_system_abstraction.jl:159-186`) decides timed vs time-free by `length(aug_state) == 3`.
  A transition from a clocked mode to a time-free one would have to consume `[x; t]` and produce
  `x`; nothing expresses that.

The *library* is closer to supporting this than the wrapper is: `_mode_physical_and_clock`
(`hybrid_symbolic_builder.jl:149-157`) already returns `nothing` for an unclocked mode, and
`HybridSymbolicModel` already handles per-mode clocked/unclocked coordinates
(`hybrid_symbolic_model.jl:56-66`). The blocking piece is transition arity.

**Fix.** P6 — make the clock a per-mode property and the reset map a per-transition
(source-arity → target-arity) object. Note that the biped and the thermostat are both
clock-free, so this too is a generalization without a present consumer.

### F7 — Dense pair tables are sized by the *global* alphabet **[read]**

The discrete solvers allocate `falses(nstates, nsymbols)` (`discrete_systems.jl:60-71`). For a
hybrid model `nsymbols = Σₖ mₖ + T` — every mode's inputs plus every switch — while any given
state can only ever use its own mode's `mₖ` columns plus the switches leaving it.

With `M` modes of comparable input count the table is ≈ `M×` larger than the information in it.
A 4-mode model with 100 inputs each and 10⁶ states allocates ~50 MB of `BitMatrix` to hold ~13 MB
of signal. `sparse_input = true` exists for the reachability and LTL solvers but **not** for the
safety solver, which is the one the hybrid safety path always uses.

**Fix.** P4 — a block-structured pair table, or route the hybrid path through `sparse_input` and
add that option to the safety solver.

### F8 — The clock lift materializes `ntime` copies of every transition **[read]**

`lift(::ClockLift, base)` (`clock_lift.jl:47-68`) replicates each base transition once per time
step. A mode with 10⁴ cells, 10 inputs and a 50 s horizon at `time_step = 0.1` (`ntime = 501`)
produces ~5×10⁶ states, and — at a branching factor of ~3 successors per (cell, input) pair —
~1.5×10⁸ transitions **for that mode alone**. This is the practical ceiling on timed hybrid models
today, and it is a product that is never inspected as a product. Arithmetic, not measured.

**Fix.** P4 — a lazy clock-lifted model whose `pre`/`post` compute the time index arithmetically
instead of storing it. The `AbstractSymbolicModel` interface is already the right seam.

### F9 — Smaller items **[read]**

| | What | Where |
| :--- | :--- | :--- |
| a | The last clock slice has no outgoing transitions (loop is `1:(ntime-1)`), so running out of time is an unmodelled deadlock rather than a stated deadline. Before `91f336f0d` those states were wrongly kept in the invariant set; they are now correctly removed. The consequence is not a tightening of the specification but the removal of an over-claim — a timed safety problem now treats the clock bound as the deadline it always was. Worth documenting, since nothing in the DSL says so. | `clock_lift.jl:55` |
| b | Modes are abstracted sequentially, though they are perfectly independent and each already has a parallel backend available. | `hybrid_symbolic_builder.jl:84-103` |
| c | A missing `time_step` on a clocked mode reaches `ClockAbstraction(sys, nothing)` and raises a `MethodError` instead of naming the mode. | `hybrid_symbolic_builder.jl:152` |
| d | `collect(tmin:tstep:tmax)` silently truncates the time domain when the span is not a whole multiple of the step. | `clock_abstraction.jl:36` |
| e | Dropped switches (empty guard intersection, reset image outside the target mode) are `@warn` only; a model can lose all of its switches and still report `OPTIMAL` for whatever the disconnected fragment admits. | `hybrid_transition_assembly.jl:70-124` |
| f | `@assert !isempty(transition_list)` gives an assertion failure rather than a diagnosis when nothing was built. | `hybrid_transition_assembly.jl:169` |
| g | An empty abstract target or safe set is not distinguished from a genuinely infeasible one — both surface as `LOCALLY_INFEASIBLE`. Hit this while writing the probes for this audit. | `optimal_control_problem.jl:141` |
| h | `_mode_time_active` compares `[1.0;;] == A` by exact float equality, and only for 1×1. | `hybrid_system_abstraction.jl:194` |
| i | `state_cost` and the continuous-time horizon are passed to the abstract problem untranslated, with `TODO`s. A horizon on a hybrid model is converted with the *model-level* `time_step`, which per-mode steps contradict. | `optimal_control_problem.jl:152-154` |
| j | JuMP cannot build a hybrid abstraction alone — `build_hybrid_problem` errors without a spec — though the MOI entry accepts `AlternatingSimulationProblem`. | `lower_hybrid.jl:401` |

---

## 4. The plan

Two rules set the order, and both were violated by the first draft of this document:

1. **A phase must not need a decision that has not been made.** A phase carrying an open API
   question cannot be the one you start on Monday.
2. **A phase must not depend on a bug a later phase fixes.** Reach-and-stay takes a safe set, and
   F3 currently makes every unspecified mode unsafe — so shipping reach-and-stay before the
   specification semantics are fixed delivers a feature crippled by a known defect.

The result front-loads everything with a present consumer and gates the rest on one appearing.

### P0 — Refuse to lose a specification silently ✅

Zero decisions, zero behaviour change for any model that is currently correct.

- `build_hybrid_problem` errors, naming the kind and the mode, on any spec kind it does not handle
  (F1). Note that a `supports_problem` declaration would *not* have caught F1: the problem never
  becomes a `ReachAndStayProblem`, so the solver's own rejection at
  `hybrid_helpers.jl:74` never fires. The check has to live in the lowering.
- Documentation truth for the switching contract (F2): say that a switch is an input the synthesis
  chooses. **Docs only** — no API change yet, see decision 1.
- Name the mode in the missing-`time_step` error (F9c); replace the empty-transition-list
  assertion with a diagnosis (F9f); document the clock deadline (F9a).

*Touches:* `lower_hybrid.jl`, `modes.jl` (docstring), `README.md:320`,
`hybrid_symbolic_builder.jl`, `hybrid_transition_assembly.jl`. *Size:* small. *Risk:* none.

### P1 — Specification semantics ✅ (except F5)

Per-kind defaults for a mode with no spec — empty for `FINAL`, the mode's whole state set for
`ALWAYS` (F3). Initial *sets* rather than the centre point, keeping the simulation start point
separate (F5). Distinguish an empty abstract target/safe set from genuine infeasibility (F9g).

**Done.** F3 (`_hybrid_spec` fills an absent mode with its own state set, for `ALWAYS` only),
F9g (`_check_nonempty` in all three hybrid sub-solvers), and F5 (`_hybrid_initial_spec` builds a
mode-indexed `HybridSpec`, discretized `OUTER` so a `start = v` point still picks up its cell).

F5 turned up a compatibility case the audit had missed: the biped hand-builds its problem with a
raw `(x, mode)` augmented point rather than through the front-end, and that is a legitimate form
— it is what `PR.satisfies` and the closed loop take. `_abstract_initial_states` accepts both,
so the front-end gains sets without breaking direct-MOI callers.

*Touches:* `lower_hybrid.jl` (`_hybrid_spec`, `_hybrid_initial_state`), `specifications.jl`,
both hybrid `build_abstract_problem`s, `results.jl`. *Size:* small–medium.
*Risk:* F3 changes results for models that wrote `Always` on a subset of modes — see decision 2.
**Open detail:** making an unspecified mode safe removes the only current way, accidental as it
is, to say "never enter this mode". `Always(LazySets.EmptySet(n))` is the obvious replacement but
has not been checked through `states_satisfying`; verify before committing to the default.

### P2 — Reach-and-stay on hybrid systems ✅

The issue's headline gap. Every piece exists: `OPDS.OptimizerReachAndStayProblem` runs in
`O(E + nm)` since `0ba8e2d95`, `states_satisfying` already evaluates a `HybridSpec` for both the
target and the safe set, and `stay_on_first_entry` already survives parsing
(`parse_scoped.jl:233`).

*Touches:* a new `OptimizerReachAndStayProblem` in `HybridSystemAbstraction/` (~60 lines, mostly
the same field/`reset!`/`MOI.optimize!` boilerplate as `safety_problem.jl`), one
`control_solver_for` method, one `build_abstract_problem`, one line in `solver_selection.jl:45`,
one branch in `build_hybrid_problem` reading `EVENTUALLY_ALWAYS`. *Size:* small — call it ~100
lines plus tests, most of it boilerplate rather than design. *Risk:* low. *Depends on:* P1.

### P3 — Report, or refuse, what was dropped ✅

Dropped switches, guard cells with no `INNER` representative, resets leaving the target domain,
resets that are not lattice-exact (F9e). **Design question I ducked in the first draft:** a report
is the weaker option. A model that silently loses *every* switch and then answers `OPTIMAL` for a
disconnected fragment should arguably fail, not warn. Proposal: hard-fail when a transition loses
all of its switches, report the partial losses, and make both inspectable as a raw attribute so a
test can assert on them.

*Size:* small. *Risk:* the hard-fail may break a model that is currently limping. That is the
point, but it belongs in release notes.

### P4 — Scale ✅ (measured, then F9b only)

Lazy clock lift (F8); block-structured or sparse pair tables, plus `sparse_input` for the safety
solver (F7); parallel mode abstraction (F9b), which is nearly free given the existing backends.

**Measured first, as this phase required.** A ring of `M` two-dimensional modes, `0.25` grid:

| modes | build s | synth s | nstates | nsymbols | pair table | pairs used |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 6.24\* | 0.46 | 450 | 52 | 0.00 MB | 37.5 % |
| 4 | 1.55 | 0.00 | 900 | 104 | 0.01 MB | 18.7 % |
| 8 | 3.25 | 0.00 | 1800 | 208 | 0.05 MB | 9.2 % |

\*includes first-call compilation.

F7's arithmetic is **confirmed exactly** — the usable fraction of the pair table halves as the
mode count doubles, tracking the predicted `1/M` — but the absolute size is 0.05 MB at eight
modes, while building the modes dominates the wall clock and synthesis rounds to zero. So F7 and
F8 are real and neither binds at any scale reachable today; **F9b is what was worth doing**, and
it is what was done. Revisit F7/F8 when a model appears whose pair table is measured in
gigabytes — the arithmetic in §3 says roughly 10⁷ states across eight modes.

**Done:** F9b, opt-in via `parallel_modes`. Opt-in rather than default because a mode's own
optimizer may use one of `Symbolic`'s threaded backends, and nesting the two oversubscribes.
Validated by asserting an identical abstraction — states, inputs and the full transition set —
against the sequential path, sharing included.

Promoted above the modelling generalizations because these three have real consumers today — the
biped is the largest hybrid model in the repository — and because all three are
behaviour-preserving, so each can be validated by asserting an identical winning set against the
current implementation, exactly as the reach-and-stay complexity fix was. **Measure first:** the
F7/F8 numbers in §3 are arithmetic on the data structures, and nothing here should be optimized
before a profile confirms which of the two actually binds.

*Size:* medium each, mutually independent. *Risk:* low.

### P5 — Sound reset maps

F4. **First task is a counterexample**, not a fix: a model with an off-lattice reset whose
synthesized controller demonstrably fails on the concrete system. Without it there is no evidence
this is a defect rather than a conservatively-worded warning, and no test that can confirm the
fix. Then over-approximate the reset image of the whole cell, reusing `System.approximation/`,
keeping the lattice-exact fast path and making the current behaviour an explicitly-named opt-in.

*Size:* medium. *Risk:* medium — transition counts grow for non-exact resets, which is the correct
cost of soundness but will surprise anyone benchmarking.

### P6 — Per-mode clocks

F6: the clock becomes a property of a mode, the reset map a per-transition map from source arity to
target arity. Multiple independent clocks are a further step and should not share this phase.

*Size:* large. *Risk:* medium-high. **No present consumer** — the thermostat and the biped are both
clock-free. Gated on decision 3.

### P7 — Forced switching, and mode invariants

F2's second half plus F2b. A transition gains a semantics — *may* (today), *must/urgent*,
*environment-chosen* — and the guard's inclusion mode flips to `OUTER` for the forced directions,
since a cell that *might* intersect the guard must carry the switch. F2b is the harder half: an
urgent semantics needs a mode invariant, which does not exist, so this phase has to introduce one
before it can define "must switch before leaving".

*Size:* large. *Risk:* high — needs its own soundness argument and a worked example with an
independently known answer. **No present consumer, and one active counter-example**: the biped
explicitly wants may-semantics. Gated on decision 3.

### P8 — Co-safe LTL over modes

Product of the hybrid automaton with the specification DBA, labels over augmented states.
Deliberately last: it is the only item needing machinery that does not exist in some form, and it
is worth little until P0–P2 make the simpler specifications trustworthy.

---

## 5. Decisions needed

1. **When to rename the switching value.** P0 fixes only the documentation. The API change —
   accepting `ControlledSwitching`, and erroring on `AutonomousSwitching` with "forced switching is
   not implemented" — should land *with* P7, not before it, so the API churns once rather than
   twice. If you would rather flip now and let existing models error with a migration message,
   that is a one-line change to P0.
2. **Is F3 a fix or a breaking change?** Making an unspecified mode fully safe changes results for
   any model that wrote `Always` on a subset of modes. Every hybrid test in the repository writes
   it on *all* modes (`hybrid.jl:397-398`, and the loop at `:193`), so nothing in-tree breaks —
   which is also why the bug survived. My reading is that any such model is currently getting an
   answer its author did not intend.
3. **Are P6 and P7 wanted at all?** This is the question the first draft asked and the repository
   already answers. The biped's own docstring (`walking.jl:72-74`) says the switch is *deliberately*
   offered to the synthesis rather than forced, and neither real model uses a clock. So there is no
   driving example for either phase, and one active argument against P7's premise. Recommendation:
   **stop after P4**, and treat P5–P7 as gated on a model that needs them. That leaves P0–P4 as the
   whole useful programme, and it is a much smaller programme than the issue implies.

---

## 6. Test plan

The existing hybrid coverage (`test/wrapper/hybrid.jl`, 11 testsets;
`test/optim/HybridSystemAbstraction/`, 3 files) tests that things *work*. It has no test that
asserts something is *rejected* for the reasons above, which is why F1 and F3 survived — checked,
not assumed: no hybrid test mentions `EventuallyAlways`, and every `Always` in the suite is
written on all modes (`hybrid.jl:193`, `:397-398`). Each phase adds:

- **P0:** an unhandled spec kind on a mode errors by name. The regression test for F1 is the
  lowering comparison used in this audit — the problem type must now *differ* with and without the
  `EventuallyAlways` constraint, where today it does not.
- **P1:** a two-mode model with `Always` on one mode only — the other mode's states are safe. An
  initial set spanning several cells — feasibility answers for the set, not for its centre.
- **P2:** the same `EventuallyAlways` model now lowers to a `ReachAndStayProblem` and solves, with
  `stay_on_first_entry` distinguishable in the winning set.
- **P3:** a model whose guard discretizes to nothing fails, rather than reporting `OPTIMAL` for the
  disconnected fragment.
- **P4:** identical winning sets before and after, on the thermostat and the 4-D biped.
- **P5:** the counterexample comes first — an off-lattice reset whose current controller fails on
  the concrete system. If that model cannot be built, F4 is not a defect and the phase is dropped.
- **P7:** a forced switch the controller would decline if it could, with a hand-computed winning
  set. This is the phase that most needs ground truth rather than a self-consistency check.

---

## 7. Not audited

The hybrid test suite was **not run** as part of this audit — the findings come from lowering
models and inspecting the resulting automata, not from a full-suite baseline. Establish one before
P0.

`PCLFBisimulationQuotient` (the other hybrid solver family), the RigidBodyDynamics extension, and
the biped problem library beyond confirming that it exercises this solver. The plotting and
dashboard path was read only where it touches `channelled_trajectory`.
