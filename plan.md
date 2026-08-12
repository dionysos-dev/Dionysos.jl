# Trajectory generation + certification pipeline — status, analysis, and plan

Working plan for: integrating the results of Florentin's master thesis (PR #569) into
the mainline pipeline, fixing the algorithmic problems found in both codebases,
restructuring the certifier into a clean two-direction family (**backward** and the new
**forward** variant, §4.4), and delivering three convincing demos.

The PR is a source of *ideas and measured evidence*, not a target to imitate. Proposals
marked **★** go beyond both codebases. Grounding: full read of the mainline pipeline
(`src/optim/trajectory_generators/`, `src/optim/trajectory_certifiers/`,
`src/system/transition_synthesis.jl`, `src/system/affine_approximation.jl`,
`ext/DionysosSymbolicsExt.jl`, the four research drivers, the tests) and a full survey
of PR #569 (396 files; ~2.3 kLOC library code, ~21.8 kLOC benchmark drivers, binaries).

---

## 0. Goals, non-goals, acceptance targets

**Goals.**
1. A **modular** pipeline: generators, certifiers, costs, providers, backends are
   swappable behind small interfaces; backward/forward certification are two thin
   instantiations of shared parts; every soundness gate is a reusable, unit-tested
   function.
2. An **efficient** pipeline: every stage measured, the known hot spots removed
   (§4.5); wall-time targets below.
3. A **sound** pipeline: no step is accepted on solver status strings or unchecked
   linearization boxes; the certificate connects to the problem's initial/target sets.
4. **Three demos, few files**: simple pendulum swing-up, double pendulum swing-up,
   articulated vehicle with obstacles — one shared driver + three configs + one
   Literate example. Not the PR's twelve 400–1200-line drivers.
5. **Evidence-driven**: every "benchmarks decide" claim in this plan has a named
   campaign in `research/` (§6b); robustness numbers are multi-seed success rates,
   never single runs.

**Non-goals (now)**: MPC/receding-horizon deployment, GPU sampling, LTL specs
(`trajectory_success(::CoSafeLTLProblem)` is a hard `false` — loud TODO), covariance
adaptation research.

**Acceptance targets** (measured by the P0 bench harness, single machine, Clarabel):
- Pendulum end-to-end (seed + MPPI + certification, formal gates on) **≥ 5×** faster
  than the P0 baseline; target < 15 s.
  *Status 2026-08-11: the qualitative bar is cleared — the baseline never certified
  at all; the pipeline now fully certifies the swing-up (§8.1) in ~79 s end-to-end.
  The wall-time target stands for the remaining optimization passes.*
- Vehicle forward maneuver end-to-end < 5 min with all formal gates passing.
- MPPI sampling throughput (rollouts/s on the pendulum) ≥ 10× baseline single-thread,
  ~linear thread scaling.
- Fast test suite stays fast: pipeline unit tests (kernels, gates, generators on the
  integrator) in the fast tier; only end-to-end demos tagged `:slow`.

---

## 1. The idea

```
problem ──▶ seed generator ──▶ MPPI refinement ──▶ nominal (x̄₀..x̄_K, ū₀..ū_{K-1})
             (abstraction /       (sampling)                 │
              RRT / NLP)                                     ▼
                                             ellipsoidal funnel chain (bwd or fwd)
                                             one small SDP per transition
                                                             │
                                                             ▼
                                        (E_1, κ_1, …, E_K, κ_K, E_{K+1})  +  guarantee:
                                        x ∈ E_k, u = κ_k(x)  ⟹  x⁺ ∈ E_{k+1}
```

- The **generator** is cheap, global, unsound (MPPI = sampling-based stochastic
  optimization; a seed supplies the homotopy class).
- The **certifier** is expensive, local, sound: per transition it linearizes the
  discrete-time dynamics, bounds the Taylor remainder with interval Hessian bounds,
  and solves one SDP (Corollary 1 of arXiv:2204.00315) synthesizing an affine feedback
  and an ellipsoid, robust to remainder and disturbance vertices.
- The product is a **funnel**: certified tube + time-varying feedback — a
  correct-by-construction tracking controller. Discrete-time, LMI-based cousin of
  funnel/LQR-tree synthesis (Majumdar & Tedrake) and DIRTREL.

---

## 2. The two codebases

### 2.1 PR #569 (Florentin)

- His library code (`heuristic_generator/`, `symbolic_certifier/`,
  `certified_solver.jl`) is an earlier, less clean version of the architecture now in
  mainline. **The mainline port is structurally better** (interfaces + MOI driver, LMI
  centralized in `ST.transition_synthesis`, LazySets conventions, tests).
- His benchmark folder is where the value is: working parameters for 4 systems, honest
  ablation reports (§9), and ideas the port did not keep (§2.3).
- Discard: 12 copy-pasted drivers (70–90 % shared), 179 PDFs + 31 GIFs + 20 MP4s +
  19 PNGs in git, accented dir names (break Windows checkouts), `println("2")` left in
  `uniform_grid_abstraction.jl`, the `get_control` shim, a double-pendulum
  `jacobian_bound` he self-flags as unreliable.

### 2.2 Mainline (`jc_issue_8`)

| Component | File | State |
| :--- | :--- | :--- |
| Generator interface | `trajectory_generators/trajectory_generators.jl` | ok |
| Abstraction seed generator | `optimizer_trajectory_generator.jl` | ok |
| MPPI generator | `mppi_trajectory_generator.jl` | works; gaps (§3) |
| Composite (seed→refine) | `composite_trajectory_generator.jl` | hardcoded to MPPI (§5) |
| Certifier interface | `trajectory_certifiers/trajectory_certifiers.jl` | ok |
| Ellipsoidal backward certifier | `ellipsoidal_backward_trajectory_certifier.jl` | 1034-line monolith; soundness gaps (§4.2) |
| Uniform-grid tube certifier | `uniform_grid_trajectory_certifier.jl` | ok — not in the PR; genuinely different certificate |
| LMI kernels | `src/system/transition_synthesis.jl` | verified correct (§4.1); caveats §4.2 |
| Affine approx + Hessian bounds | `affine_approximation.jl`, `ext/DionysosSymbolicsExt.jl` | sound; **recompiles per call** (§4.5) |
| MOI driver | `trajectory_certification_optimizer.jl` | one-shot, no feedback (§6) |
| Experiments / tests | `research/TrajectoryCertificationOptimizer/`, `test/optim/` | Integrator + simple pendulum only |

### 2.3 What the PR has that mainline lost in the port

1. **Periodic unwrap before certification** (+ `lift_state_near_reference` replay).
   Mainline just avoids the seam (`research/.../simple_pendulum.jl:323-326`). Must port.
2. **Terminal set = inscribed ellipsoid of the target box + containment check.**
3. **Certifier-aware MPPI cost** (terminal-ellipsoid distance in cost/success/truncate).
4. **Input reserve** (plan on shrunk `U_plan`, certify on full `U`). Mainline's
   pendulum script got the direction backwards (§4.2-E).
5. **Statistical validation harness** (`run_kappa_statistical_check`).
6. **Reach-avoid audit** (funnels vs `X` and obstacles). Mainline never checks.
7. **`σ̃ = tan δ`** vehicle input reparametrization (control-affine; still unported).
8. **`symbolic_system(...)` constructors in `problems/`** (only `NonLinear` has one).
9. **Working parameters + measured ablations** — §8, §9.

---

## 3. MPPI — errors, missing features, better optimizer

References: Williams et al. 2016/2017/2018; SMPPI; log-MPPI; MPPI-Generic. What we
have is *iterated open-loop MPPI with elitism* — the α = 1 corner of IT-MPC (no
importance-sampling correction) used as a one-shot trajectory optimizer. Same in the PR.

### 3.1 Genuine errors

- **E1 — Temperature/cost-scale mismatch (likely the main performance killer).**
  Research drivers use `λ = 1.0` against costs of scale O(10²–10³): `exp(−ΔS/λ)` puts
  all weight on the argmin sample — the update degenerates to random search, invisible
  because nothing logs the effective sample size. **Fix**: ESS-adaptive λ (choose λ per
  iteration so `ESS = 1/Σw²` hits a target fraction; bisection, a few lines); log ESS
  + weight entropy.
- **E2 — Truncation makes crashing cheap** (`mppi_trajectory_generator.jl:218-227`).
  On violation the rollout breaks: violating state invisible to the cost, shorter
  rollout accumulates less cost; empty-input `sum` throws on first-step violation.
  Latent only because every call site passes `hard_constraint = false` against a
  `true` default. **Fix**: never truncate inside sampling — finite penalty, keep
  rolling (or freeze); full-horizon comparable costs.
- **E3 — Projection bias**: rollouts use `proj(u+ε)`, the update averages raw `ε`.
  **Fix**: average `ε_eff = u_roll − u_nom`. One line.

### 3.2 Missing (ports & standard practice)

- **M1 — Structured noise + optional IS term**: `GaussianMPPINoise(Σ)` (closure kept
  as fallback); with it, `γ·Σ_t ū_tᵀΣ⁻¹ε_t`, `γ = λ(1−α)`, default `α = 1` documented.
- **M2 — Smoothness pressure**: `Δu`/`Δ²u` penalties (certifier-friendly); optional
  Savitzky–Golay on the update.
- **M3 — Cost builders**: `TrackingCost`, `TerminalEllipsoidCost`, `InputEffort`,
  `InputSmoothness`, `DomainPenalty` — composable structs, not 100-line closures.
- **M4 — Input reserve**: project onto shrunk `U_plan`, certify on true `U`.

### 3.3 ★ Beyond the PR

- **★S1 — Stage-cost rollout kernel.** Both codebases allocate a `Trajectory` per
  sample and re-collect references per sample per iteration. Replace the hot path with
  `stage_cost(x,u,k)` + `terminal_cost(x)` evaluated in a preallocated, type-stable
  kernel (no `Vector{Any}`, no per-sample allocation). Expected ~10× before threads.
- **★S2 — Threaded + antithetic sampling**: `Threads.@threads` with per-thread buffers
  and split RNG streams (deterministic given seed); mirrored ±ε pairs ≈ halve the
  samples needed.
- **★S3 — CEM update mode**: for one-shot optimization an elite-refit update often
  beats softmin and adapts σ for free; same infrastructure —
  `update_rule = :mppi | :cem`, benchmarks decide.
- **★S4 — σ annealing; spline-parametrized noise** (knots → smooth candidates, lower
  search dimension). Stretch.
- **★S5 — Abstraction value function as terminal cost**: Dionysos-native — the coarse
  abstraction's cost-to-go as `φ(x_H)` lets the horizon shrink well below task
  duration and guides globally through obstacles. Stretch, but the biggest strategic
  lever for the double pendulum.

---

## 4. Certification — soundness, the forward variant, efficiency

### 4.1 Verified correct

I re-derived the Schur complements of all S-procedure blocks in
`transition_synthesis.jl` — correct. The remainder model is sound: per-output interval
Hessian eigenvalue bound over the linearization box (`DionysosSymbolicsExt.jl:67-96`);
`Lip·(δx + δu)` with `δx ≥ ‖Δx‖²`, `δu ≥ ‖Δu‖²` is the 2nd-order Taylor bound with the
½ dropped (2× conservative). Cosmetics: lipschitz vector allocated `nx+nu+nw` long,
only `1:nx` filled (per-*output* — rename/resize); Δw-dependent remainder unmodeled
(assert `ΔW = 0` until it is).

### 4.2 Soundness & correctness gaps (ranked)

- **A — Fixed mode has no box-consistency gate (unsound).** Required radii are logged,
  never enforced; the SDP cap is `maxδx` (30 / 10⁶ in scripts), not the linearization
  box (±0.2). Florentin measured the consequence (margins 0.044/0.029 on a
  "successful" run). **Fix**: mandatory gate (`:inconsistent_box ⟹ :infeasible`), plus
  optional in-LMI cap `δx ≤ (min_i δx_box_i)²` (sound: source centered at `x_k`).
- **B — Scaled Lipschitz transform under-estimates the remainder.** `|T⁻¹|·L` maps the
  error vector but leaves `(δx+δu)` in z-units; truly `‖Δx‖² ≤ σ_max(T)²·δx_z`. The
  PR's vehicle scalings (σ_max ≈ 15) under-estimate up to **225×** — matching his
  finding that no vehicle config passed formal checks at box scale 1.0. **Fix**:
  compute bounds on the scaled dynamics — free under global normalization (§4.3).
- **C — Funnels never constrained to `X` / away from obstacles**: the certificate is
  not reach-avoid. Post-hoc gate now (port PR audit); in-LMI half-spaces where the
  parametrization allows (backward: SOC in `L`; forward `:free`: SOC in `Q₂`, §4.4).
- **D — Certificate not connected to the spec**: no `E_terminal ⊆ target_set` check,
  no initial-set coverage report. Port inscribed-ellipsoid builder + containment gate;
  report coverage margins. (Forward closes the initial side by construction.)
- **E — Pendulum experiment certifies against `U_cert = ±10.5` on a ±7 plant** (and
  MPPI can emit excluded deadband inputs). Correct direction: `U_plan ⊆ U`,
  `U_cert = U`. Use the PR's convex `benchmark_up_convex` objective for the showcase.
- **F — Cost matrix passed where a factor is expected**: LMIs bound `‖Λ·[x;u;1]‖²`;
  callers pass the PSD matrix — bounds `vᵀS²v`, an *under*-estimate when `S ≺ I`.
  Today's call sites are accidentally idempotent. **Fix**: factor internally
  (cholesky); non-idempotent regression test.
- **G — Funnel collapse has no floor**: PR double pendulum dies at k=35 with radii
  3·10⁻⁷ — legal, meaningless, then fatal. **Fix**: optional `L ⪰ r_min·I`; see also
  the maximin objective (§4.6).
- **H — Periodic states**: unwrap-before-certify + lifted replay, default when the
  problem declares periodic dims.
- **I — Default traps**: `EllipsoidalBackwardOptions.λ = 1.0` deletes the volume
  objective (kernel default 0.01; every experiment overrides); MPPI
  `hard_constraint = true` (every call site overrides). Defaults must match real use.
- **★K — Solver status is not a certificate.** Backward accepts
  `ALMOST_OPTIMAL / NEARLY_FEASIBLE_POINT`; an approximately-feasible point is not
  sound. **Fix**: `ST.validate_transition` — plug the returned solution into the PSD
  blocks, check `eigmin ≥ 0`; if marginal, deterministically shrink the source
  (`L ← (1−η)L`, bisect η) and report the certified shrink. Cheap dense eigchecks;
  makes soundness backend- and tolerance-independent. Lives in `ST` beside the kernels
  (it needs the block algebra); the certifier's gate just calls it.
- **J — Small stuff**: `skip_already_in`/`skip_point` dead flags; forward kernel
  demands `OPTIMAL` while backward accepts `ALMOST_OPTIMAL` (harmonize via ★K);
  `K = F/Lval` on near-degenerate `L` (regularize or fail loudly); degenerate noise
  boxes enumerate 2^nw duplicate vertices; `_collect_kappas` returns `Nothing[]`;
  `set_problem!` runs twice per solve in the MOI driver.

### 4.3 Normalization — yes it helps; do it globally

PR ablations (pendulum): **no scaling ⟹ fail at k=52/66; with scaling all 66 certify**;
`trajectory_std` won his sweep. Mechanisms: SDP conditioning; the volume objective
otherwise buys metres and starves radians; the LMI's `δx` is a Euclidean budget that
anisotropic units misallocate. But the per-step transform is where bug B lives, and
mainline's `state_scaling` is dead code (no call site sets it).

**Design: normalize once, globally.** Nondimensionalize state/input at problem setup
(`x̃ = T⁻¹x`, `T` = characteristic scales) and build the symbolic discrete-time
dynamics in normalized coordinates. Bug B disappears by construction; Hessian bounds
are computed in working coordinates; boxes / `maxδx` / `r_min` / trace objectives get
uniform meaning; the abstraction seed conditions better too. Per-step scaling survives
only as an advanced option with the `σ_max(T)²` correction.

### 4.4 The forward certifier variant — core feature, not an extra

Direction, formally:
- **Backward (current)**: target known, *source shape is a decision variable* → the
  `L`-congruence trick, a volume objective, adaptive boxes — all because the remainder
  depends on a box the SDP has not chosen yet.
- **Forward**: entry ellipsoid `E_1` given (smallest shape-fixed ellipsoid covering
  `problem.initial_set`, centered at `x̄_1` — closes gap D by construction); for
  `k = 1..K` synthesize `κ_k` and the target `E_{k+1}`.

**Why the LMIs get simpler:**

1. **Source is data ⟹ the linearization error is known before solving.** `δx` is a
   constant (`λ_max(Q_k)`), the Hessian bound is computed over the exact box of `E_k`:
   the x-side consistency gate holds **by construction** — no adaptive boxes, no
   growth/safety/scale search. Only the u-side involves the unknown `κ_k`; the
   existing input-proximity block bounds it with one `δu` variable (or a fixed budget
   + one posterior check).
2. **The target can be free — even its shape — and the SDP stays convex.** The target
   enters the reach block only as `Q₂ = P₂⁻¹`, alone in the (3,3) position, never
   multiplying `K`: linear. This is the classical asymmetry of ellipsoid calculus — a
   free *target* is convex as-is; a free *source* is what forced backward's machinery.
   Two modes sharing everything:
   - `target_shape = :fixed(P̂)` — scale `α` free (`Q₂ = α·P̂⁻¹`), objective `min α`.
     Cheapest; yields the contraction profile.
   - ★ `target_shape = :free` — `Q₂` a decision variable, objective `min tr(Q₂)` or
     `min λ_max(Q₂)`. **Caveat**: true volume `min logdet(Q₂)` is
     concave-minimization — not representable; trace/λmax are the convex size
     measures and, in normalized coordinates, the right ones. This mode is
     one-step-optimal ellipsoidal tube propagation with feedback synthesis: shapes
     adapt to the local dynamics. Thin tubes are *not* a pathology here (tighter
     certificate ⟹ smaller remainder next step); the only hazard is numerical —
     step k+1 inverts `Q₂` — so add a conditioning sandwich
     `q_min·I ⪯ Q₂ ⪯ q_max·I`.
   - Bonus of the `Q₂` parametrization: domain/obstacle half-spaces become
     `aᵀQ₂a ≤ (b − aᵀc₂)²` — SOC in `Q₂` — **reach-avoid enforced inside the LMI**.
     Free center `c₂` possible (enters affinely) but off by default (re-anchors the
     next linearization).
3. **`α` is a per-step contraction certificate** (fixed-shape mode); the chain output
   becomes a contraction profile with an `α ≤ α_max` fail-fast guard. Free-shape mode
   reports `tr(Q₂)` per step. Default `:free` for production, `:fixed` for
   diagnostics and the duality test.

**When to use which** (both stay; direction is a certifier option):

| | Backward | Forward |
| :--- | :--- | :--- |
| Answers | largest certified *entry* set for a given target | given the entry set to serve, where the tube ends |
| Spec connection | terminal exact; entry coverage reported after | entry exact; `E_{K+1} ⊆ target_set` checked at the end |
| Failure | mid-chain infeasibility | late but informative (see where the tube inflates) |
| Machinery | adaptive boxes, volume objective | none of it |
| Re-planning | re-plan the **prefix** into the certified suffix | re-plan the **suffix** from the certified prefix |

- ★ **Bidirectional handoff**: forward pass from the initial set + backward pass from
  the target; certified as soon as the forward tube at some `k` is contained in the
  backward funnel at `k` (one `UT.is_included` per step). Splits long chains, reuses
  both partial results. Cheap once both directions exist — belongs in the loop (§6).

Kernel work: `solve_transition_forward` (both modes) beside the two existing kernels,
from the same shared S-procedure blocks — the shared-block design was made for this.

### 4.5 Efficiency budget — where the time goes, what removes it

**P0 baseline (measured 2026-08-10, Clarabel, single thread, fixed rng —
`bench/trajectory_pipeline.jl`; `bench/results/` is gitignored so the reference lives
here):**

| System | seed | mppi | provider cold/warm | chain | status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| integrator | 4.0 s | 0.39 s (2555 rollouts/s) | 11.6 s / 1.3 ms per call | 21.7 s (1671 ms/step, incl. JuMP+Clarabel JIT) | certified |
| pendulum | 8.1 s | 0.86 s (1157 rollouts/s) | 2.4 s / 5.5 ms per call | 140 ms/step | **failed_k = 19** |

Two P0 lessons that correct earlier assumptions: (i) RuntimeGeneratedFunctions cache
compiled code, so warm provider calls cost only the per-call symbolic
re-differentiation (1.3–5.5 ms) — the precompilation lever is a 2–4× on the pendulum
chain, not 10×; the multi-second "cold" numbers are one-time JIT. (ii) The pendulum
certification **fails at k=19 deterministically** under the baseline config — a
reproducible instance of the reported "pipeline not performing well", now the standing
test case for P1–P5.

PR measurements (single-thread Mosek): LMI chain dominates everywhere; the abstraction
seed is 5–50× the MPPI cost.

| Stage | PR measured | Lever | Expected |
| :--- | :--- | :--- | :--- |
| Seed (abstraction) | 1.9 s (pendulum) … 161 s (vehicle) | deliberately **coarse seed grid**; RRT seed (`UT.search/RRT` exists); ★ NLP collocation generator | 10–100× |
| MPPI sampling | 1.6 … 30 s | ★S1 stage-cost kernel; ★S2 threads + antithetic; E1 ESS-λ (fewer iterations wasted) | ~10× + thread scaling |
| Provider (inside chain) | hidden in LMI time | **★ stop recompiling**: `build_affine_approximation` re-derives `Symbolics.jacobian`, every Hessian entry, and re-runs `build_function` **per call** (`DionysosSymbolicsExt.jl:258-292`); the adaptive loop calls it up to `max_iters + |scales|` times per step. Compile once at provider construction; per call = interval evaluation | plausibly the single biggest chain win; measure in P0 |
| SDP chain | 28.6 s (pendulum) … 316 s (vehicle) | **forward mode** (no adaptive iterations at all); `:first_consistent` default backward (his ablation: certifies everything the 5× line-search does); ★ ball remainder model (`‖e‖ ≤ ‖Lip‖·(δx+δu)`, one multiplier instead of 2^nx vertex blocks — default at nx ≥ 4); parallel box candidates (backward); maximin objective (no logdet cone) | 2–20× |
| Model build | unmeasured | profile JuMP build vs solve; direct MOI only if build dominates | measure first |

Backend stays injected. **Decision: Clarabel only for now** (a Mosek license exists but
is deferred — the PR's Mosek settings stay documented for the day it is enabled). This
raises the value of two items: the maximin objective (no logdet cone — Clarabel reaches
`LogDetConeTriangle` only through exponential-cone bridging) and the ball remainder
model (Clarabel feels SDP size more than Mosek). All §0 targets are Clarabel numbers.

### 4.6 ★ Remaining structural ideas — challenged

- **★ Maximin-radius objective (backward)**: logdet/trace volume is gameable by
  pancake ellipsoids — the observed collapse mode. `max r s.t. L ⪰ r·I` (min
  semi-axis) prevents collapse by construction, needs no exotic cone, and turns the
  §4.2-G floor into the objective. Default for backward; logdet kept as option.
- **TVLQR-shaped funnels — demoted.** Originally pitched for speed (scale-only SDP)
  and shape continuity. The forward `:free` mode now delivers adapted shapes from the
  SDP itself, and forward already removed the adaptive-box cost. TVLQR shapes remain
  interesting only if *backward* stays the bottleneck (they'd give it a cheap `P̂`
  and a scale-only variant). Revisit after P0/P3 numbers; do not build on spec.
- **★ Two-pass remainder tightening** (re-solve once with the remainder recomputed
  from the actual synthesized shape): optional, try on the double pendulum only.

---

## 5. Architecture — modularity is a deliverable

### 5.1 Certifier family

```
src/optim/trajectory_certifiers/
  trajectory_certifiers.jl      # interface (set_problem!/set_trajectory!/certify!/…)
  ellipsoidal/
    options.jl                  # EllipsoidalChainOptions + per-direction extras
    context.jl                  # problem + trajectory + provider + normalization, built once
    steps.jl                    # backward_step! / forward_step!  → one shared StepRecord
    gates.jl                    # box-consistency, residual validation (calls ST.validate_transition),
                                # reach-avoid, terminal/initial containment — pure, direction-agnostic
    chain.jl                    # run_chain!(ctx, direction) — iteration, stop/retry, result assembly
    adaptive_boxes.jl           # backward-only, quarantined
    diagnostics.jl              # step summaries, contraction/volume profiles, out of the hot path
  uniform_grid_trajectory_certifier.jl
```

### 5.2 Generator family (same treatment — currently missing)

```
src/optim/trajectory_generators/
  trajectory_generators.jl      # interface incl. set_seed_trajectory! (promoted; error stub)
  costs.jl                      # composable stage-cost builders (M3): TrackingCost,
                                # TerminalEllipsoidCost, InputEffort/Smoothness, DomainPenalty
  rollout.jl                    # ★S1 typed rollout kernel (shared by MPPI/CEM; threads live here)
  mppi/
    noise.jl                    # GaussianMPPINoise(Σ) + closure fallback (M1)
    updates.jl                  # :mppi softmin (+ optional IS term) and :cem elite refit (★S3)
    generator.jl                # the loop: seed, iterate, elitism, diagnostics (ESS, entropy)
  optimizer_trajectory_generator.jl
  composite_trajectory_generator.jl   # generic over the interface; fallback knobs
```

### 5.3 Design rules (pipeline-wide)

1. **Variant by dispatch, not flags**: `forward_step!`/`backward_step!` are methods
   over one `StepRecord`; `run_chain!` never branches on a symbol. Same for update
   rules and noise types on the generator side.
2. **Gates are pure functions** `gate(record, ctx) → Union{Nothing, FailureReason}`,
   composed into a list, each unit-tested with an engineered failing input. Soundness
   does not live inside a 900-line stepping function.
3. **All LMI algebra stays in `ST.transition_synthesis`** — three kernels + ★
   `validate_transition`; the certifier assembles and gates, it never touches blocks.
4. **Everything expensive is built once, in constructors**: precompiled provider,
   SDP backend, normalization, rng, cost builders (references hoisted). `certify!` /
   `generate!` only evaluate.
5. **One result vocabulary**: both directions return the same
   `EllipsoidalCertificationResult` (funnel, controllers, step records, profiles);
   plotting, the statistical harness, and `FunnelController` are written once.
6. **The certified controller is a controller**: `FunnelController <: ST` dynamic
   controller (state = step index, `output_control = κ_k(x)`, `domain = E_k`), so
   `ST.get_closed_loop_trajectory`, `simulate`, plotting, and validation reuse the
   existing protocol. Demo = "one `optimize!`, then simulate like any controller".
7. **Costs are data**: builders with a `stage_cost(x,u,k)` protocol; hand closures
   remain the escape hatch, never the recommended path.

### 5.4 Seeds (modularity dividend)

The abstraction seed costs 5–50× the MPPI refinement (PR timings) for a homotopy
class. Behind the same interface: deliberately coarse seed grid; RRT seed; ★
direct-collocation generator (JuMP + Ipopt) — for smooth problems fast *and*
certifier-friendly, challenging MPPI outside obstacle-cluttered tasks. Benchmark
question, not religion.

---

## 6. The generate ⇄ certify loop

Neither codebase has feedback from certification failure into generation; the MOI
driver is one-shot (and calls `set_problem!` twice). Build in order:

1. **Retry ladder**: on `failed_k`, retry the step with grown boxes/scales — surfaced
   in the driver result, not a bare `false`.
2. **Partial-chain re-planning**: the certified part of the chain stays valid.
   Backward: re-plan the *prefix* into the entry ellipsoid `E_{failed_k+1}` (via
   `TerminalEllipsoidCost`), warm-started from the failed trajectory. Forward:
   symmetrically re-plan the *suffix* from the certified tube at `failed_k`. Iterate;
   every round reuses all certified work.
3. **Bidirectional handoff** (§4.4): run forward from the initial set and backward
   from the target, succeed on first per-step inclusion. An option of the driver once
   both directions exist — cheap and particularly promising for the double pendulum.
4. **Certifiability-aware generation (DIRTREL-flavored)**: cost terms for what makes
   steps uncertifiable — saturation margin (input reserve), `Δu`/`Δ²u`, ★ local
   nonlinearity proxy (precompiled Hessian bound at `x̄_k` × step length — free after
   the provider precompilation). No SDP in the generation loop.

---

## 6b. Evidence: benchmarks and robustness campaigns

Two tiers, two homes:

- **`bench/` — timing regression** (the P0 harness): per-stage wall times on fixed
  configs, re-run at every phase end, results committed as the reference table.
- **`research/TrajectoryCertificationOptimizer/campaigns/` — scientific campaigns**:
  robustness, ablations, and the shoot-outs this plan defers to. One **shared campaign
  runner** (grid of configs × rng seeds × systems → CSV summary + regenerated plots);
  each campaign is a ~50-line config on top of it — this is the explicit
  anti-copy-paste device the PR lacked (his 12 drivers *were* his campaigns).
  Binaries are never committed; small CSV summaries are.

**Methodology rule**: MPPI is stochastic — a configuration's result is the success
rate + cost/time quantiles over **≥ 20 rng seeds**, never one run. The PR's numbers
are single-seed; treat them as anecdotes until reproduced.

| Campaign | Phase | Question it settles |
| :--- | :--- | :--- |
| **C1 — generator ablation** | P2 | `:mppi` vs `:cem`; ESS target; antithetic on/off; success-rate vs sample-budget fronts |
| **C2 — certifier ablation** | P3 | normalization variants (reproduce Florentin's sweep on mainline: none / fixed / trajectory-std); `remainder_model = :vertices` vs `:ball` (conservatism vs time); objective maximin / logdet / trace |
| **C3 — direction shoot-out** | P4 | backward vs forward `:fixed` vs forward `:free` vs bidirectional — certified entry volume, chain time, failure position, per system |
| **C4 — loop effectiveness** | P5 | fraction of failed certifications recovered by retry / re-planning / handoff; iterations and extra time to close the chain |
| **C5 — seed comparison** | P5 | fine vs coarse abstraction vs RRT — seed time vs downstream MPPI effort and certification success |
| **C6 — disturbance robustness** | P6/stretch | certify with `W ≠ 0`; Monte-Carlo replay under process noise; `planning_input_scale` (input reserve) front |

---

## 7. Roadmap

Each phase ends green on the fast suite, formatted, one commit. Ordering principle:
**baseline first, skeleton before gates** (gates land once, in their final home — no
build-then-move churn).

**P0 — Baseline bench harness + campaign runner** *(you cannot claim "efficient" or
"robust" without them)*
- `bench/trajectory_pipeline.jl`: per-stage timings (seed / MPPI throughput / provider
  / SDP chain / gates) on the integrator and pendulum; runs in minutes; results
  committed as the reference table. Re-run at the end of every phase.
- The shared campaign runner in
  `research/TrajectoryCertificationOptimizer/campaigns/` (§6b): configs × seeds ×
  systems → CSV + plots; multi-seed success-rate reporting built in. Smoke campaign on
  the integrator to prove the runner.

**P1 — Certifier family skeleton + soundness (REF + FIX)**
- Module split per §5.1 (behavior-preserving move of the backward certifier;
  `adaptive_boxes.jl` quarantined; diagnostics extracted).
- Gates land in `gates.jl`: mandatory box-consistency (A), terminal/initial
  containment (D), reach-avoid post-hoc (C), collapse floor (G).
- `ST.validate_transition` (★K) + harmonized statuses; cost factor via cholesky (F) +
  non-idempotent test; defaults fixed, dead flags removed (I, J).
- Tests: each gate against an engineered failing input; golden values unchanged.

**P2 — MPPI overhaul**
- E1–E3 (ESS-λ, penalty rollouts, effective noise) + diagnostics.
- §5.2 layout: `rollout.jl` kernel (★S1), threads + antithetic (★S2), `noise.jl`
  (M1 + optional IS term), `updates.jl` (`:mppi`/`:cem`, ★S3), `costs.jl` (M3, M4).
- Tests: ESS regression, penalty-vs-truncation ordering, seeded determinism
  (single/multi-thread). Bench: throughput vs P0. **Campaign C1** settles
  `:mppi` vs `:cem` and the ESS/antithetic defaults.

**P3 — Certifier performance & conditioning**
- ★ Precompiled symbolic provider (differentiate/compile once at construction).
- Global nondimensionalization; per-step scaling fixed (`σ_max(T)²`) or quarantined;
  PR round-trip checks as tests (B, §4.3).
- ★ Maximin objective default for backward; ★ `remainder_model = :ball` (default
  nx ≥ 4); `:first_consistent` default; parallel box candidates. Bench: chain time.
  **Campaign C2** validates the normalization design and the new defaults.

**P4 — Forward variant**
- ★ `solve_transition_forward` (`:fixed(P̂)` α-mode and `:free` Q₂-mode with trace
  objective + conditioning sandwich + optional in-LMI half-spaces) — §4.4.
- `forward_step!` (no adaptive boxes), `α ≤ α_max` guard, contraction profile, entry
  ellipsoid from `initial_set`.
- Tests: forward golden values (integrator); the **duality property test** — on a
  linear system, backward's synthesized entry set fed to the forward certifier must
  certify with `α ≤ 1 + tol` at every step (cross-validates both kernels, both steps,
  and the shared gates at once). **Campaign C3** ranks the directions per system.

**P5 — Pipeline plumbing + loop**
- `FunnelController` + simulate/plot integration; periodic unwrap default (H);
  statistical harness as test utility.
- **Decision: the full loop ships in the first release** — retry ladder, partial-chain
  re-planning (prefix/suffix by direction), *and* the bidirectional handoff, all as
  driver modes. Composite generic over the interface + fallback knobs; coarse-grid /
  RRT seed options; driver `set_problem!` dedup. **Campaigns C4 + C5** measure loop
  recovery rates and pick the default seed.

**P6 — Problems & demos** *(few files)*
- Port to `problems/`: `σ̃ = tan δ`, `symbolic_system` constructors,
  `benchmark_up_convex` objectives.
- One shared demo driver + three configs (§8); one Literate example (pendulum); tests
  wired (`:slow` for end-to-end). Bench: final table vs P0 targets (§0).
  **Campaign C6** (disturbance robustness, input-reserve front) backs the demos'
  robustness claims with multi-seed numbers.

**P7 — Stretch**
- ★S5 abstraction-value terminal cost; ★ NLP collocation generator; in-LMI
  half-spaces for backward; ★ two-pass remainder tightening; TVLQR shapes *only if*
  backward remains the bottleneck (§4.6); spline noise / SMPPI; double-pendulum full
  chain; 2-trailer bonus.

---

## 8. Per-system battle plans (PR's proven parameters as starting points)

### 8.1 Simple pendulum — reference demo (low risk) — **ACHIEVED 2026-08-11**

`research/TrajectoryCertificationOptimizer/demo_pendulum.jl`: the driver loop
certifies the complete 48-step swing-up in round 1 (~56 s + 23 s seed, Clarabel,
fixed rng) — all soundness gates passing, terminal ⊆ target, on the plant's true
±4.5 input set (the P0 baseline failed at k≈19 against an unsound ±10.5 set),
ending in a 48-step `ST.FunnelController`. Decisive ingredients: seam-aware target
(`UT.set_in_period`), the `:cem` update (softmin averaging failed the multimodal
energy pumping), unwrap + 2π-shift via the driver's `prepare_trajectory` hook, and
the `:maximin` objective across 48 chained funnels. Gate-reported residuals:
entry-funnel coverage of the initial set (gap D — forward direction / entry
enlargement is the designed answer) and the handoff not yet nesting here.

Florentin certified 66/66 with adaptive boxes + scaling. Start: Δt=0.1, H=55–75,
`nsamples`≈3–8k (expect far fewer after P2), `niter`=20, MPPI λ=1.75 until ESS-λ, σ_u
=0.5–1.8, `planning_input_scale`≈0.75; certifier λ=0.005, terminal = inscribed
ellipsoid (shrink 0.85), normalization ≈ `[0.22, 0.25]`. Convex-X objective; fix the
`U_cert` direction (E). Run **both directions**: backward for the maximal entry set,
forward for initial-set coverage + contraction profile; the pendulum doubles as the
duality test on a real system. Deliverable: certified swing-up + funnel animation +
Monte-Carlo validation.

### 8.2 Articulated vehicle — the impressive one (medium risk)

Forward maneuver with two obstacles as headline; parallel parking second if cheap.
Needs `tan δ`, periodic handling (θ, ϕ), reach-avoid (post-hoc, or in-LMI via forward
`:free`), and P3 normalization (his scalings magnify, σ_max ≈ 15 — unsound under the
current transform; his audits only passed with 2–3× boxes). Start: Δt=0.1–0.2,
H=33–60, `nsamples`=1800 forward (~10× more reverse — antithetic + threads matter),
MPPI λ=0.3, σ=[0.45, 0.15]; certifier λ=0.001, terminal r=0.45. Success = formal pass
of reach-avoid + box + residual gates, not just chain feasibility.

### 8.3 Double pendulum — the hard one (high risk, plan the fallback)

Known failure (PR audit): funnel collapse at k≈35/66 (radii ~3·10⁻⁷), 16 error
vertices/step, huge Hessian bounds. Attack order: (a) global normalization;
(b) maximin backward / forward `:free` with conditioning sandwich; (c) ball remainder
+ precompiled provider to make the SDP count affordable; (d) shrink terminal
(john-shrink 0.5–0.7, `ΔU`≈0.02); (e) partial-chain re-planning and the
**bidirectional handoff**; (f) smoother MPPI (M2/★S4) — his σ_u=0.8, λ=1.75, 1800
samples, Δt=0.05, H=51, plus his two-stage terminal-refinement MPPI as an option.
**Fallback that is still impressive**: certify the *capture phase* (last ~20 steps
into upright) end-to-end with the abstraction handling energy pumping — an honest
partial certificate beats a fake full one.

**OUTCOME (shipped, `demo_double_pendulum.jl`)**: the capture-phase certificate, at
a much deeper level of understanding than planned. Generation solved raw-CEM (no
abstraction): resonant seed scan u = A·sin(ωt+φ) + height shaping (1 + cos θ₂) +
terminal pull — a direct ~6 s swing. Certified plant = the **RK2 map** (RK4's four
symbolic self-compositions stack-overflow the Hessian compiler — measured; RK2
compiles in ~20 s and generation ≡ certification map, zero mismatch). New machinery
shipped en route: `source_cap` SOC rows in `solve_transition_backward` +
`ChainOptions.domain_cap` (funnels confined to X and to the linearization box by
construction — the box-scale ladder becomes the funnel-size dial). The wall: a
**four-lever ablation, all null** (input headroom ×2, ω₂ pace cap, Δt 0.05→0.025 —
runway 1.35→1.43 s, Δt-invariant — and the free-angle domain: the optimal swing
never leaves θ₁ ∈ ±π/2). Funnel bleed through the ballistic ascent is per-unit-time
intrinsic (one saturated input vs four states): volumes 14.9 → 1e-12 over ~1.4 s.
Shipped: 57-step capture chain certified standalone (terminal ⊆ target, domain
gates green), reported as a depth profile — guaranteed recovery region
[0.06, 0.06, 0.3, 0.3] (rad, rad, rad/s, rad/s) at 0.12 s out, [0.014, …, 0.07] at
0.5 s out — with closed-loop FunnelController validation (50/50) from the
meaningful depth. Note vs the PR — different task: its MPPI experiment ran
`benchmark_up_convex`, the UP-UP handstand (both angles π ± 25°, torque ±8.5,
ω ∈ ±6, both angles periodic) — that is where k=34/66 collapsed; our shipped demo
is the swing_up objective. Its plant also has an operator-precedence slip
(Coriolis term not divided by l·α in both `dynamic` and the symbolic model), so
its numbers are on an internally consistent but non-textbook system, with
unenforced fixed boxes (14.8× violations per its own audit).

**UP-UP HEAD-TO-HEAD DONE (`demo_double_pendulum_upup.jl`)**: our stack on his
exact benchmark. Generation transfers with three deltas (both-links height
shaping, NO pace cost — the flip needs ω near the ±6 wall, SHORT 3 s horizon —
6 s dilutes CEM and fails, measured): a 56-step / 2.8 s double inversion, both
angles inside π ± 25° (his: 2.55 s). Certification: **24 enforced transitions**
(k=56→33, V_max ≈ 1.9, wall mid-flip at k=32 — the same per-unit-time ballistic
bleed; certified runway 1.2 s vs swing_up's 1.4 s), capture chain re-certified
standalone, depth profile ([0.03, 0.03, 0.15, 0.15] physical semi-axes at
0.25 s out), 50/50 closed-loop from the meaningful depth. Verdict vs his
32-of-66: every one of our transitions is consistency-enforced,
collapse-proofed, and domain-capped on the textbook plant — his near-terminal
(fat, valuable) ellipsoids are exactly the informal ones.
**REMAINDER-MODEL LADDER (measured, controlled — same trajectory, only the
model changes)**: on the up-up chain at Δt = 0.05, `:ball` wall k=32 (24
transitions), `:john_ball` (NEW — the box's John ellipsoid √n·diag(Lip)·B as a
shaped Petersen ball, per-axis radii at :ball's single-block cost) k=32 —
IDENTICAL, because both ball covers are corner-tight and the win lies in the
exact joint corner support — and `:vertices` k=24 (32 transitions, **+33%
runway**). On the swing_up chain at Δt = 0.025, `:vertices` leaves the runway
unchanged (wall physics-dominated at fine steps) but **doubles the guaranteed
recovery regions at every depth**. Rule: `:vertices` whenever n ≤ 4 (16
blocks/SDP is affordable); `:ball` for long chains in higher dimension
(`:john_ball` kept — sound, tested — but no measured advantage). Both DP demos
ship with `:vertices`.

**BIDIRECTIONAL VERDICT (measured, entry-scale ladder ×1.0/×0.5/×0.25/×0.1 of
I)**: the forward tube explodes at the same physical stretch where the backward
funnel bleeds — at entry ×0.1 (±0.02 rad) it survives 8 steps (~×2.7 growth)
then inflates ×7 in 3 steps through the pump ascent and dies; larger entries die
immediately (the saturated nominal leaves ~1.0 of correction headroom against
±1.7 rad/s deviations). No handoff exists at any scale. The wall is now
characterized bidirectionally with seven null ablations around it — the clean
signature of a per-step certificate on a rank-1-B plant: u defends one direction
per step, the expansion lives in the other three. Remaining research levers (not
tuning): (1) MULTI-STEP (lifted) transitions — certify E_k → E_{k+2} in one LMI
with composed maps and two gains, giving the certificate the 2-step coupling
channel u → ω₁ → ω₂ that per-step chains can never see (helps both directions);
(2) anti-greedy windowed chain optimization (~5-step joint SDPs — both audits
show per-step greed contributes); (3) adaptive Δt schedules (volume/compute win —
measured 16× terminal volume at fine Δt — but not a wall-breaker: runway is
Δt-invariant).

---

## 9. PR #569 measured results worth remembering (do not re-learn these)

- Fixed boxes "succeed" unsoundly (margins 0.044/0.029) — why gate A is mandatory.
- No-scaling fails at k=52 (pendulum); `trajectory_std` scaling won; `trajectory_range`
  fails outright.
- Vehicle: no configuration passed all formal post-hoc checks at box scale 1.0;
  `none`/`trajectory_std` became FORMAL at 2.0 — consistent with bug B.
- Monte-Carlo failures dominated by chain exits in the first ≤3 steps, not target
  misses — entry-funnel size / initial coverage (gap D), not tracking.
- Reverse vehicle maneuvers need 1–2 orders of magnitude more MPPI samples.
- LMI chain ≫ MPPI everywhere; abstraction seed 5–50× the MPPI cost.

## 10. Decisions & open questions

**Decided:**
1. **Backend: Clarabel only for now.** A Mosek license exists but is deferred; all
   benchmarks and §0 targets are Clarabel numbers (see §4.5 for the two items this
   re-prioritizes).
2. **Loop scope: everything in.** The first release ships the retry ladder,
   partial-chain re-planning (§6-2), and the bidirectional handoff (§6-3) as driver
   modes.

**Decided (was open):**
3. Double pendulum: the capture-phase certificate is the shipped demo (user-approved;
   see §8.3 OUTCOME). The full-chain wall survived a four-lever ablation and is
   characterized as intrinsic per-unit-time funnel bleed through the ballistic ascent.

**Open:**
4. Should the tube-based `UniformGridTrajectoryCertifier` get showcase treatment (a
   genuinely different certificate), or stay a tested-but-quiet alternative?
