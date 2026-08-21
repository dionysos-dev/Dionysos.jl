# Plan — Adaptive Cruise Control benchmark

> **Status — T1 is built** (2026-08-21). `problems/AdaptiveCruiseControl/`,
> `examples/AdaptiveCruiseControl/`, `test/problems/adaptive_cruise_control.jl`, and the
> growth-bound testset all land and pass. T2–T4 are untouched, by decision. Two findings from
> building it are recorded in §11; the first is a soundness bug this benchmark exposed in the
> safety solver, now fixed.

Add **Adaptive Cruise Control (ACC)** to the Dionysos benchmark library: a longitudinal
car-following problem where an ego vehicle tracks a desired cruising speed but must never
violate a time-headway constraint with respect to the vehicle in front.

It is worth adding for three reasons that no current benchmark covers at once:

1. It is the **canonical correct-by-construction case study** in the symbolic-control
   literature (Nilsson et al., IEEE TCST 2016), so results are comparable with published work.
2. Its safe set has a **closed-form maximal controlled-invariant set** (§6), so the abstraction
   can be checked against ground truth rather than against itself — none of the existing
   `problems/` entries offer that.
3. Its natural formulation carries an **adversarial exogenous input** (what the lead vehicle
   does), which is the one modelling ingredient the library currently has no example of. §5
   shows it is expressible today, but not by the route one would first reach for.

---

## 1. Answers to the three framing questions

| Question | Answer |
| :--- | :--- |
| **What kind of system?** | **Continuous-time, nonlinear, single-mode.** Nonlinearity is the quadratic aerodynamic drag. It is *not* natively hybrid; a genuine hybrid variant exists (§4, T3) but it is a different modelling choice, not the same model written differently. |
| **What specification?** | **Safety is the core** (never violate the headway), but the honest full ACC specification is **reach-and-stay under a safety constraint**: `□ safe ∧ ◇□ cruise`. Both are directly expressible — `PR.SafetyProblem` and `PR.ReachAndStayProblem` (which carries `target_set` *and* `safe_set`). A third, reach-avoid, gives a time-to-cruise value function. |
| **State / input dimension?** | **2 states, 1 input** in the reference formulation `x = (z, v)`, `u = a`. **3 states, 1 input** once the lead's speed is a state rather than a parameter. **5 states, 1 input** for the actuator-lag (ARCH-COMP) version, which is beyond the uniform grid. |

**Recommendation: start at 2-D.** It is the published benchmark, it is cheap enough to iterate
on in seconds, and it is the version with analytical ground truth.

---

## 2. The model

Ego vehicle as a longitudinal point mass with rolling resistance and drag; `z` is the bumper-to-
bumper gap to the lead vehicle, `v` the ego speed, `v₀` the lead speed.

```
ż = v₀ − v
v̇ = ( −F_r(v) + u ) / m ,        F_r(v) = c₀ + c₁ v + c₂ v²
```

Working in **acceleration units** (`a = u/m`) instead of wheel force keeps every quantity on a
readable scale and makes the input grid step interpretable as m/s²; do that.

Jacobian and growth bound (hand-written; `L[i,j] ≥ |J[i,j]| off-diagonal, L[i,i] ≥ J[i,i] on it`):

```
J = [ 0   −1                    ]        L = [ 0    1                       ]
    [ 0   −(c₁ + 2 c₂ v)/m      ]            [ 0   −(c₁ + 2 c₂ v_lo)/m     ]
```

`L` is constant in `u`, and `L₂₂` is evaluated at the **lowest** speed in the state box, where
`J₂₂` is least negative. `L₁₂ = 1` says the gap uncertainty grows by `τ ·` the speed uncertainty
— which is the whole coupling in this system.

### Parameters and where each number comes from

Following the practice established for the biped document: state the provenance, so the reviewer
knows which numbers are literature and which are ours to argue about.

| Symbol | Value | Provenance |
| :--- | :--- | :--- |
| `m` | 1650 kg | Ames et al. 2017, benchmark value |
| `c₀, c₁, c₂` | 0.1 N, 5 N·s/m, 0.25 N·s²/m | Ames et al. 2017 |
| `v₀` (lead speed) | 13.89 m/s (50 km/h) | Ames et al. 2017 |
| `v_d` (desired) | 24 m/s (86 km/h) | Ames et al. 2017 — deliberately **above** `v₀`, which is what makes the constraint bind |
| `a_max` | 0.3·g = 2.943 m/s² | Ames et al. 2017 (`c_a = c_d = 0.3`, a wheel-friction bound) |
| `τ_h` (time headway) | 1.8 s | Ames et al. 2017; also the ISO 15622 default — **chosen**, and the single parameter the answer is most sensitive to |
| `X` state box | `z ∈ [0, 100] m`, `v ∈ [0, 30] m/s` | **ours.** §6 shows the invariant set needs `z ≲ 74 m` at `v = 30`, so 100 m contains it with margin |
| `τ` (period) | 0.1–0.5 s | **ours**, to be tuned in M2. ACC control periods are typically 50–100 ms; the abstraction may want longer |
| `h_z, h_v` | 0.5 m, 0.1 m/s (start) | **ours**, to be tuned in M2 |
| `h_a` | 0.2 m/s² (≈ 30 inputs) | **ours** |

A note that will matter when tuning: with `h_a = 0.2` and `τ = 0.5`, `τ·h_a = 0.1 = h_v`, so the
speed axis is an **exact lattice** — each input translates the speed grid onto itself. The gap
axis is not (`τ·h_v = 0.05 ≠ h_z`), and the drag term breaks exactness anyway. Worth checking
whether aligning both axes is affordable; it would make the abstraction deterministic.

---

## 3. The specifications, mapped onto `PR.ProblemType`

| # | Informal | Dionysos type | What the certificate is |
| :--- | :--- | :--- | :--- |
| S1 | Never tailgate: `□ (z ≥ τ_h v)` | `PR.SafetyProblem` | maximal controlled-invariant set — comparable to the closed form in §6 |
| S2 | Reach cruising speed and hold it, never tailgating | `PR.ReachAndStayProblem(sys, I, T, S)` | winning set + memoryless controller |
| S3 | From a slow, distant start, reach the cruise band without tailgating | `PR.OptimalControlProblem(...; safe_set)` | value function = worst-case time to cruise |

For **S2**, the target set is a union, and the union is the point: the real ACC objective is
*"cruise at `v_d`, **or** follow the lead safely if `v_d` is unreachable"*.

```
T = { |v − v_d| ≤ ε }  ∪  { |v − v₀| ≤ ε  ∧  z ≥ τ_h v + margin }
```

Built with `UT.set_union`; `MP.get_states_from_set` accepts any bounded `LazySet`.

> **Expect S2 with the single-component target `{|v − v_d| ≤ ε}` to be infeasible**, because
> `v_d = 24 > v₀ = 13.89` — you cannot cruise faster than the car in front forever. That is not a
> bug to work around; it is the benchmark proving something true, and the two-component target
> is the fix. Both runs belong in the example page: the infeasible one is the more instructive.

---

## 4. Variants, in the order they should be built

| Tier | State | Input | Time | What it exercises | Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **T1** | 2: `(z, v)` | 1: `a` | continuous | `GROWTH` + hand-written Jacobian bound; S1/S2/S3 | **build first** — cheap, published, has ground truth |
| **T2** | 3: `(z, v, v_l)` | 1: `a` | continuous | bounded-disturbance lead ⇒ genuinely nondeterministic abstraction (§5) | build second — this is the scientifically novel part |
| **T3** | 2 + mode | 1 | hybrid | `@mode` / `add_transition!`; ego switches speed-control ↔ gap-control (ISO 15622 architecture) | build third, if T1/T2 land cleanly |
| **T4** | 5: `(z, v, γ, v_l, γ_l)` | 1 | continuous | first-order actuator lag; the ARCH-COMP NNCS model | **document as out of scope** for the uniform grid; note it as a lazy/multi-scale target |

T4 is listed so the omission is deliberate and recorded, not an oversight. `x_lead` and `x_ego`
collapse into `z`, which is why the nominal 6 states reduce to 5.

---

## 5. The adversarial lead needs care — two traps

This is the part of the design most likely to produce a silently unsound benchmark, so it is
written out rather than discovered during implementation. Both traps were checked against the
current source.

**Trap 1 — the lead's acceleration is not an input.** It is adversarial: the synthesis must
quantify *universally* over it. Putting it in the input vector makes it controller-chosen, and
the resulting "guarantee" would be worthless.

**Trap 2 — hybrid modes are also controller-chosen.** In the Dionysos front-end the mode is
selected by the controller (see the hybrid section of the DC-DC example, where the switch *is*
the control). So "the lead switches between cruise / brake / accelerate modes" is **not**
expressible as a hybrid model with adversarial modes. This is exactly why T3 is framed as the
*ego* switching control regimes, which genuinely is a controller decision.

**What does work, with no change to the toolbox.** The synthesis already quantifies universally
over `Post(q, u)`, so *any* nondeterminism in the transition relation is treated adversarially.
It is therefore enough to widen the reachable set to cover every admissible lead behaviour. Note
that `ContinuousTimeGrowthBound(system; jacobian_bound = L)` integrates `ṙ = L(u)·r` — purely
multiplicative, with **no additive term** — so passing a `jacobian_bound` cannot express this.
Two routes that can:

- **Preferred — pass `growthbound_map` directly** instead of `jacobian_bound`. The attribute is
  already accepted by the optimizer and short-circuits the derivation, and its contract is
  `(r, u, tstep) -> r'`. Integrate `ṙ = L(u)·r + w̄` with `w̄ = (0, 0, ā_l)`, reusing
  `ST.runge_kutta4`. This is the standard SCOTS-style formulation and stays inside `GROWTH`.
- **Fallback — `approx_mode = USER_DEFINED`** with an `overapproximation_map`, if the reachable
  set ever needs a shape the radius formulation cannot express.

**M4 must include an explicit soundness test for this**, not just a growth-bound check: Monte-
Carlo the true reachable set under randomly sampled lead-acceleration signals and assert
containment in the abstraction's `Post`. A wrong `w̄` fails silently in exactly the way the
existing `test/problems/growth_bounds.jl` was written to catch for `L`.

---

## 6. Analytical ground truth (T1) — the validation criterion

With the lead at constant `v₀`, the barrier `h = z − τ_h v` obeys `ḣ = (v₀ − v) − τ_h a`. Under
maximum braking `a = −a_max`, `ḣ = v* − v` with

```
v* = v₀ + τ_h · a_max              (= 13.89 + 1.8 · 2.943 = 19.19 m/s)
```

So `h` can only decrease while `v > v*`, and integrating over the braking transient down to `v*`
gives the maximal controlled-invariant set in closed form:

```
z ≥ τ_h · v                                       for v ≤ v*
z ≥ τ_h · v + (v − v*)² / (2 a_max)               for v > v*
```

At `v = 24` this needs `z ≥ 47.1 m`; at `v = 30`, `z ≥ 73.9 m`. The boundary is piecewise
linear-then-parabolic in the `(v, z)` plane, which also makes for a figure worth plotting.

Two caveats to state on the page rather than bury:

- This **ignores drag**, which assists braking, so it is an **inner** approximation — the true
  set is slightly larger (~5 % in `a` at `v = 24`, where `F_r/m ≈ 0.16` m/s²).
- It is a **continuous-time** result. Sampling at period `τ` can only shrink the winning set.

Hence the acceptance criterion for M2: the abstraction's invariant set **contains** the closed
form up to one cell, and is contained in it inflated by the sampling margin. A computed set that
*exceeds* the analytical one by more than the drag correction means the abstraction is unsound.

---

## 7. Files to add and touch

Following the checklist in `CLAUDE.md` §5 and `docs/src/developers/examples.md`.

**New**

| Path | Contents |
| :--- | :--- |
| `problems/AdaptiveCruiseControl/adaptive_cruise_control.jl` | module `AdaptiveCruiseControl`: `Params`, `dynamic`, `jacobian`, `jacobian_bound`, `growthbound_map` (T2), `system`, `safety_problem`, `reach_and_stay_problem`, `optimal_control_problem`, `invariant_set_closed_form`, `system_plot!`. Mirror `problems/DCDC/dcdc_converter.jl` exactly. |
| `problems/AdaptiveCruiseControl/adaptive_cruise_control_hybrid.jl` | T3 only, once T1/T2 land. Mirrors the `Thermostat/` multi-file layout. |
| `examples/AdaptiveCruiseControl/adaptive_cruise_control.jl` | runnable MOI driver, DCDC-style: build the abstraction once, then solve S1, S2, S3 against the cached abstraction. |
| `test/problems/adaptive_cruise_control.jl` | see M2/M4 acceptance criteria. |
| `docs/src/examples/jump/adaptive_cruise_control.jl` | Literate page through the JuMP front-end. |

**Touch**

| Path | Change |
| :--- | :--- |
| `test/problems/growth_bounds.jl` | add an `@testset "AdaptiveCruiseControl"` — mandatory, this file covers every hand-written bound in `problems/` |
| `test/runtests.jl` | add `("./problems/adaptive_cruise_control.jl",)` to `TEST_FILES` (untagged if it stays under a few seconds) |
| `docs/make.jl` | add `"adaptive_cruise_control"` to `ORDER["jump"]`; without it the page still builds but is appended alphabetically |
| `docs/src/references.bib` | §10 |
| `bench/definitions.jl` | optional, M6: add an ACC case to the cross-solver suite |

`examples/README.md` needs no edit — it describes the folder convention, not a list.

---

## 8. Milestones

Each milestone ends with `julia -e 'using JuliaFormatter; format(".")'` and one commit, per the
repo convention (`[ACTION] module: description`, lowercase, ≤ 60 chars).

**M1 — the model.** `problems/AdaptiveCruiseControl/adaptive_cruise_control.jl` with T1 dynamics,
Jacobian, growth bound, `system`, and the three problem constructors. Add the growth-bound
testset.
*Done when:* `julia --project=test test/problems/growth_bounds.jl` passes with the new testset.
*Commit:* `ADD problems/acc: adaptive cruise control model`

**M2 — safety, validated against ground truth.** Solve S1 through `UniformGridAbstraction`.
Implement `invariant_set_closed_form` and the containment test of §6. Tune `τ`, `h_z`, `h_v`,
`h_a` here — record what was tried and what was kept.
*Done when:* the computed invariant set brackets the closed form as specified, and the whole
solve runs in under ~30 s.
*Commit:* `ADD problems/acc: safety problem and invariant set`

**M3 — reach-and-stay and reach-avoid.** S2 with both the single-component (infeasible) and
two-component (feasible) targets, and S3. Closed-loop trajectories for each.
*Done when:* S2-single reports infeasible, S2-union reports feasible, and every simulated
trajectory satisfies `trajectory_success` for its own problem.
*Commit:* `ADD problems/acc: reach-and-stay and reach-avoid specs`

**M4 — the adversarial lead (T2).** Third state `v_l`, custom `growthbound_map`, plus the Monte-
Carlo containment test of §5.
*Done when:* the containment test passes and the invariant set is visibly *smaller* than T1's —
if it is not, the disturbance is not reaching the abstraction and something is wrong.
*Commit:* `ADD problems/acc: bounded-disturbance lead vehicle`

**M5 — driver and documentation.** `examples/…` driver, `docs/src/examples/jump/…` Literate page
with the animated dashboard (`system_plot!` closure + `animate_trajectory_dashboard`, per the
house convention), `ORDER` entry, bib entries.
*Done when:* `julia --project=docs/ docs/make.jl` is green. Watch for the `@ref` /
`CurrentModule` trap — a page citing a function by `@ref` resolves against `Main` unless the page
declares `CurrentModule`; the `@cite` keys must exist in `references.bib`.
*Commit:* `ADD docs/acc: adaptive cruise control example`

**M6 (optional) — T3 hybrid and the bench case.**
*Commit:* `ADD problems/acc: hybrid ego-mode variant`

Run `julia --project -e 'using Pkg; Pkg.test(; test_args = ["--fast"])'` at the end of M4 and the
full suite before the PR. Run Julia from PowerShell, not Git Bash.

---

## 9. Open decisions

These change the work materially and are worth settling before M2 rather than during it.

1. **Is the headline benchmark T1 or T2?** T1 is the published, checkable one; T2 is the one that
   demonstrates something the library cannot currently show off. The plan builds T1 first either
   way, but which one the documentation page leads with is a positioning choice.
2. **Sampling period `τ`.** A realistic ACC period (50–100 ms) may make the abstraction nearly
   self-looping at a readable grid resolution; a longer `τ` (0.5 s) buys a cleaner abstraction and
   an exact lattice on the speed axis but is no longer a realistic controller period. Pick one and
   say which, in the provenance table.
3. **Force or acceleration as the input.** Acceleration is recommended (readable grid, directly
   comparable to `a_max`), but force is what the source paper uses, so the numbers in a
   side-by-side comparison would need converting.
4. **Does `τ_h = 1.8 s` stay fixed, or become a swept parameter?** A sweep over `τ_h` producing a
   family of invariant sets would be a strong figure, and cheap at 2-D.

---

## 10. References to add to `docs/src/references.bib`

- **`nilsson2016correct`** — Nilsson, Hussien, Balkan, Chen, Ames, Grizzle, Ozay, Peng, Tabuada,
  *"Correct-by-Construction Adaptive Cruise Control: Two Approaches"*, IEEE Transactions on
  Control Systems Technology 24(4):1294–1307, 2016. The direct precedent: ACC solved by
  abstraction-based synthesis, with the lead as a bounded disturbance.
- **`ames2017control`** — Ames, Xu, Grizzle, Tabuada, *"Control Barrier Function Based Quadratic
  Programs for Safety Critical Systems"*, IEEE Transactions on Automatic Control
  62(8):3861–3876, 2017. Source of the model and every parameter in §2.
- **`iso15622`** — ISO 15622:2018, *Intelligent transport systems — Adaptive cruise control
  systems*. For the time-headway convention. Optional; cite only if `τ_h` is discussed as a
  standards-driven choice rather than a benchmark constant.

Verify each entry against the actual publication before committing — a wrong bib entry passes CI
silently.

---

## 11. What building T1 turned up

### 11.1 A soundness bug in the safety fixed point — found, fixed

`compute_largest_invariant_set` kept states with **no outgoing transition at all**. Such a state
is one every input drives out of the domain, so it cannot be part of any controlled-invariant
set. The counter-based fixed point removes a state only when a pair being *disabled* drops its
count to zero, and a count that starts at zero never drops — so the state survived, its
predecessors survived with it, and the controller reported itself *defined* there while having
no control to offer.

Minimal reproducer: three states, one input, `1 → 1`, `2 → 1`, and `3` with no transition. The
invariant set came back `{1, 2, 3}` instead of `{1, 2}`.

On this benchmark it was not subtle. Every cell near the top of the gap axis at low speed is
trapped — the gap opens at up to 13.6 m/s and no input closes it before the state leaves the
box — so the invariant set ran to `z = 98` at *every* speed. After the fix the upper boundary
rises with speed, as it must: `[2, 34]` at 0.25 m/s, `[20, 90]` at 10 m/s, `[42, 98]` at 20 m/s.

Fixed in `src/optim/discrete_systems/safety_problem.jl` and, identically, in
`_invariant_with_floor` in `reach_and_stay_problem.jl` (the same shape; floor cells stay
protected, and `_compute_seed` already refuses to win a state on an empty `post`, so they always
have a pair). Regression test in `test/optim/discrete_systems/safety.jl` — the existing test
there walks `domain(contr)` and so could never have seen it.

The reachability solvers are **not** affected: their counter enables a pair only on a decrement
event, and a pair with no transitions is never in any `pre` list, so it is never enabled.

### 11.2 Reach-and-stay was `O(depth ⋅ E)` — fixed, ~125× here

The benchmark's first measurements had the union-target reach-and-stay taking **379 s**
(`stay_on_first_entry`) and **437 s** (`◇□`) on the driver's abstraction — against 6.6 s for
reach-avoid on the same automaton. Neither `stay_on_first_entry` nor `early_stop` moved it much,
which pointed at the algorithm rather than the configuration.

Both solvers grew the winning set one `Pre` layer per round, and computed each layer by
re-walking every transition: `O(depth ⋅ E)`. Depth is a property of the **plant**, not of the
specification — near the lead's speed this benchmark's gap closes by centimetres per step, so
its reachability runs hundreds of layers deep and every one of them cost a full sweep. The
reach-avoid solver never had this problem: its counter-based attractor touches each transition
once overall.

Both now use that same device (`_attract!` in `reach_and_stay_problem.jl`). Measured, same
abstraction, no `early_stop`:

| | before | after |
| :--- | ---: | ---: |
| reach-and-stay, `stay_on_first_entry` | 379 s | **3.0 s** |
| reach-and-stay, `◇□` | 437 s | **3.5 s** |

For `stay_on_first_entry` this is a strict complexity improvement to `O(E + nm)`: the invariant
core is computed once, so the approach is a single attractor. For `◇□` the outer μ/ν loop
remains, but the round count drops from the reachability depth to the number of times the
invariant core actually grows, and the attractor work is amortised across all rounds by keeping
one counter table.

**Verification.** A/B against the previous implementation, both injected into the same process
so they run on the identical automaton: on the pendulum swing-up (33 480 states, 4 864 670
transitions, periodic and genuinely nondeterministic) both variants return **identical** winning
sets and identical sets of controlled cells, 40× and 14× faster. Every other reach-and-stay user
in the repo re-run and green (§11.3).

That A/B earned its keep: the first version of the rewrite was **unsound**, returning 30 716
winning states where the reference returned 29 914. A cell won by the attractor in one round can
join the invariant core in the next, and it was then re-seeded as a frontier target, decrementing
its predecessors' counters a second time — dropping some pair's count to zero while a successor
was still unwon. The fix separates "already relaxed" from "winning"; a counter may be decremented
once per transition, ever. `test/optim/discrete_systems/reach_and_stay.jl` pins it with an
automaton that returns `[1, 2, 4]` under the bug and `[2, 4]` correctly.

### 11.3 Reach-and-stay coverage

That solver had **no discrete-level test file** — only the abstraction-level
`test/optim/UniformGridAbstraction/reach_and_stay.jl`. Added
`test/optim/discrete_systems/reach_and_stay.jl`: hand-built automata where the winning set is
worked out by hand, covering both readings of "and stay" and where they diverge, the safe set,
targets with no outgoing transition, nondeterminism, the double-relaxation regression, and
`early_stop`.

Every reach-and-stay user in the repo, re-run after the change (times include ~20–35 s of Julia
startup and JIT, and were taken under CPU contention):

| | |
| :--- | ---: |
| `test/problem/problems.jl` | 19.4 s |
| `test/wrapper/specifications.jl` | 72.0 s |
| `test/optim/UniformGridAbstraction/reach_and_stay.jl` | 22.6 s |
| `docs/src/examples/getting_started.jl` | 136.7 s |
| `docs/src/examples/jump/integrator_ltl.jl` | 301.2 s |
| `examples/Pendulum/simple_pendulum_reach_and_stay.jl` | 66.7 s |
| `examples/Thermostat/thermostat_system.jl` | 38.4 s |

A warning for anyone re-measuring these on Windows: PowerShell's `$?` is `$false` for a native
command that merely wrote to stderr, so `julia ... *> $null` reports a clean run as a failure.
Use `$LASTEXITCODE`.
