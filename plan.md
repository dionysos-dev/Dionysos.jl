# Plan — Integration of Riccardo's velocity-control biped work

> **Status (2026-08-17): implemented, all phases.** `problems/BipedRobot/4D_model/robot_problem.jl`
> (model, sound carving, IK target, commensurable grids), `src/mapping/grid_mapping/lattice.jl`
> (`intersample_safe_time_step`, `is_lattice_exact`), the measured inter-sample jump `@warn` in the
> UGA solver, `src/optim/discrete_systems/bounded_input_variation.jl` (pair-keyed Dijkstra +
> `BoundedInputVariation` attribute, lifted concrete→symbols by the UGA front-end), the laptop driver
> `examples/BipedRobot/biped_4d_velocity.jl`, and tests (`test/problems/biped_robot_4d.jl`,
> `test/mapping/lattice.jl`, `test/optim/discrete_systems/bounded_input_variation.jl`).
>
> Findings from the implementation worth keeping:
> - **Sound carving proves Riccardo's resolution was too coarse**: at `dx = 0.1` the Lipschitz margin
>   (~5.5 cm) provably disconnects the free space around his obstacle — the "obstacle jumps" he saw
>   were the unsound center-test masking an infeasible discretization. At `dx = 0.05` (margin
>   ~2.7 cm, `|θ| ≤ 1.2`) the footstep over an 8 cm × 3 cm step is feasible: BFS-connected, and the
>   synthesized controller crosses in 52 steps (the BFS optimum).
> - **Memory bound, not time bound**: 3⁴ = 81 inputs × 2.5 M cells ≈ 200 M transitions exceed laptop
>   RAM; restricting to one-joint-per-step inputs (`state_input_filter`, 9 effective inputs, 21 M
>   transitions) keeps the whole run under 30 s of abstraction time on 14 threads.
> - The slew-rate-limited synthesis reaches the same 52 steps with the measured max input variation
>   exactly at the bound — the constraint is active and satisfied, with rest-to-rest ramps.

Branches analyzed: `robot-velocity-control`, `state-constraints` (superset of the former),
`robot-graph-algorithm`, `robot-abstraction-test` (already absorbed into master — nothing to take).
All are based on pre-refactor-v0.2 master (June 2026), so everything must be ported to the
current API. Scope decided with Julien: **laptop simulation only, no cluster work** (the generic
SLURM infra already exists under `problems/BipedRobot/execution/`).

---

## 1. What Riccardo did

### 1.1 The velocity-control model (`robot-velocity-control` → `state-constraints`)

`problems/BipedRobot/4D_model_vcontrol/` (~1100 lines):

- **4-D kinematic model**: states = 4 joint angles (left/right hip and knee), inputs = joint
  angular velocities, dynamics `ẋ = u`. The torque level is delegated to the motors' low-level
  velocity controllers. This halves the dimension vs the 8-D torque model (and removes the
  stiff multibody integration entirely — `f(x,u)` is trivial instead of a 0.1 s RigidBodyDynamics
  roll-out per transition).
- `geometry.jl`: forward/inverse kinematics of the two-leg mechanism (stance/swing legs,
  grounded-foot convention).
- **Physical obstacle in the domain model** (`state-constraints`): a Cartesian obstacle (a
  triangular step, `LazySets.Polygon`, in the swing-foot plane) is pulled back into joint space
  by `remove_infeasible_cells` (cell-by-cell preimage test), together with the ground
  constraint (swing foot `y ≥ 0`). The state domain becomes `X_full ∖ ⋃ infeasible cells`.
- `compute_target_set`: inverse kinematics sweep producing the set of joint configurations
  that place the swing foot at the next footstep location (a 2-D manifold gridded into cells).
- Two-phase solver use = the canonical abstraction-caching pattern
  (`AlternatingSimulationProblem` to build the abstraction once, then swap in the
  `OptimalControlProblem`), with `ThreadedBackend`, JLD2 save/load, and a gif animation.

**The hidden gem (implicit in his constants, never stated):** he chose the input grid
*commensurable* with the state grid — `du = 10·dx`, `τ = 0.1` ⟹ `τ·du = dx`. Every step
translates the state grid exactly onto itself, so the abstraction is **exact and deterministic**
(a bisimulation, not just an alternating simulation): zero spurious transitions. The
discretization-nondeterminism curse — the main killer for the 8-D torque abstraction — vanishes
by construction. This is the scientific core of the velocity-control approach and must be made
explicit, checked, and documented in the integration.

### 1.2 Input-rate-constrained synthesis (`robot-graph-algorithm`)

A backward Dijkstra in `src/optim/discrete_systems/` enforcing `d(u⁻, u) ≤ Δ` between
consecutive inputs along the path, plus a final target input `target_u`. Since inputs are
velocities, this is an **acceleration limit**: it makes the synthesized velocity profile
trackable by the motors' low-level controllers — exactly the assumption the whole approach
rests on. Wired in by adding 3 fields to `OptimalControlProblem` and a branch in the discrete
OCP solver. Includes a 3×3 toy example (`low_acceleration_path.jl`).

---

## 2. Assessment — is it nice? Is it the right direction?

**Yes on both counts.** The direction is sound and well-matched to Dionysos:

- Velocity-level control is the standard answer to abstraction-based control of
  high-dimensional mechanical systems: the abstraction lives where the curse of dimensionality
  is manageable (4-D kinematics) and the certified/uncertified boundary is explicit (the
  low-level velocity tracking is the uncertified layer — its tracking error can later be
  re-imported as a disturbance bound, see §5).
- The commensurable-lattice design makes the abstraction *exact* — the strongest possible
  soundness relation, at the cheapest possible cost (`f` is a translation). Very few
  abstraction papers get to work with an actual bisimulation.
- The obstacle work is the right formulation: physical constraint in Cartesian task space,
  pulled back to the state-space domain model, solved with unchanged solver machinery.
- The input-rate Dijkstra closes the loop with the physical motors. The three pieces together
  form a coherent story: *reduce → abstract exactly → synthesize with actuator-aware
  constraints*.

**But the code cannot be merged as-is.** It predates the v0.2 refactor, and it has one
correctness gap (the "obstacle jump" Riccardo observed himself) plus several bugs — all
fixable, see §3.

---

## 3. Issues found (and what the fixes are)

### 3.1 The obstacle jump — root cause and the correct sampling-time bound

Audit result: the abstraction pipeline computes, per (cell, input), **only the reach set at
time τ** (`transition_kernels.jl:43`); nothing anywhere constrains the trajectory during
`[0, τ]`. On top of that, two problems specific to Riccardo's code:

1. **Center-only infeasibility test**: a cell is removed iff the *center's* swing-foot position
   hits the obstacle. The obstacle's joint-space preimage is a thin shell (≈ 5 cm obstacle /
   ≈ 0.37 m kinematic Lipschitz ≈ 0.13 rad ≈ one cell at `dx = 0.1`), so it crosses cells
   without containing their centers → those cells survive → the abstract trajectory passes
   "through" the obstacle. Unsound.
2. `CENTER_SIMULATION` is an **under-approximation** mode (the `src/` docstring misleadingly
   says "very conservative"). In the exact-lattice design it happens to be exact, but nothing
   verifies the lattice property.

**Which approx mode to use? → `CENTER_SIMULATION`, guarded by an explicit lattice check.**
When the lattice property holds (`τ·u ∈ hx·ℤ` for every input-grid point), simulating the
center is *exact*: the cell's true reach set is exactly one cell, and the center lands in it —
the abstraction is a bisimulation, deterministic, one transition per (cell, input). `GROWTH`
with `jacobian_bound = 0` is **not** an acceptable substitute here: the reach box then
coincides exactly with a cell, and the OUTER index-range discretization uses closed boundaries
(`get_pos_lims_outer`, `tol = 0.0`, `grid.jl:59-70`) — a face-aligned box picks up all its
neighbors, i.e. **3ⁿ = 81 successors instead of 1**. Sound, but it destroys both determinism
(killing the input-rate Dijkstra) and tightness (worst-case over 81 successors). So:
`CENTER_SIMULATION` + an `assert`-style commensurability check at problem construction
(`τ·du ≈ hx` within fp tolerance, input grid centered on multiples of `du`) is *both* the fast
and the sound choice — the check is precisely what upgrades it from unsound-in-general to
exact-here. If a user breaks commensurability, fail loudly with a message pointing to `GROWTH`.

**Riccardo's suggested fix (τ ≤ obstacle_width / v_max) is the right instinct but not sound**:
it prevents jumping fully *over* the obstacle, yet still allows *grazing* — between two free
endpoints the continuous path can clip an obstacle corner. The correct bound is
grid-referenced, not obstacle-referenced:

> **No-jump theorem.** If (i) cell removal is **OUTER** — every cell *intersecting* the
> obstacle preimage is removed — and (ii) **τ · u_max,i ≤ hx_i on every axis** (at most one
> cell of displacement per axis per step), then the continuous segment between consecutive
> states stays inside `source_cell ∪ target_cell`, both certified obstacle-free. No
> intermediate-time check is needed.
>
> Hence the sampling-time upper bound: **τ_max = min_i hx_i / u_max,i.**
> The exact-lattice design satisfies it with equality.

Exact OUTER preimage removal is hard (nonlinear kinematics), but a Lipschitz
over-approximation suffices and stays cheap: remove the cell when
`dist(foot(center), obstacle) ≤ L_FK · half_diagonal`, with `L_FK` the segment-length-weighted
Lipschitz constant of the forward kinematics. (Precedent in-repo: the tube certifier's
anti-gap rule `max_step ≤ 0.5·rmin`, `uniform_grid_trajectory_certifier.jl:185`.)

### 3.2 Point bugs in his code

| Where | Bug | Fix |
| :--- | :--- | :--- |
| `remove_infeasible_cells` | `HyperRectangle(rec.lb*0.95, rec.ub*0.95)` scales toward the **origin**, not the cell center — removed boxes are shifted for any cell away from 0 | shrink around the cell center (or drop the shrink entirely once removal is OUTER) |
| `remove_infeasible_cells` | O(n⁴) sweep over ~10⁶–10⁷ cells | decompose: swing foot = hip(θ₁,θ₂) + swing leg(θ₃,θ₄) → two O(n²) precomputations |
| `compute_target_set` | same cell pushed hundreds of times (1e-3 IK sweep, no dedup); `inflated_rec` computed then dead (`rec` pushed instead) | accumulate a `Set` of cell indices; decide inflation deliberately |
| Dijkstra (`compute_optimal_controller_bounded_var`) | PQ keyed on `(q, u)` but value table and controller keyed on `q` alone — paths reaching `q` with different `u` overwrite each other → sub-optimal and possibly incomplete | true DP on the product `(q, u_prev)` (§4.3) |
| Problem struct | 3 untyped fields bolted onto `OptimalControlProblem`; `target_u === Nothing` compares to a *type* | leave `OptimalControlProblem` untouched; constraint becomes a solver attribute |
| Simulation script | `CENTER_SIMULATION` used with no justification — unsound in general, exact here only by an unchecked accident of the constants | keep `CENTER_SIMULATION`, add the explicit lattice-exactness assertion (see §3.1) — do **not** switch to `GROWTH` (3ⁿ face-aligned fanout) |

### 3.3 API port (pre-refactor → current)

| Old (his branches) | Current |
| :--- | :--- |
| `UT.HyperRectangle(lb, ub)` | `LazySets.Hyperrectangle(; low, high)` (keyword! positional = center/radius) |
| `UT.LazySetMinus` / `UT.LazySetUnion` | `UT.set_minus` / `UT.set_union` (grid layer consumes them natively, hole mode inverted → sound) |
| `UT.get_dims` | `LazySets.dim` |
| `x_traj, u_traj = get_closed_loop_trajectory(...)` | single `ST.Trajectory`; `ST.states(traj)` / `ST.inputs(traj)` |
| `"efficient"` attribute | removed — setting it now throws |
| grid/mapping/problem/backend APIs | unchanged (`GridFree`, `get_states_from_set`, `AlternatingSimulationProblem`, `ThreadedBackend`, all UGA attributes) |

---

## 4. Integration plan

### Phase 1 — `problems/BipedRobot/4D_model/robot_problem.jl` (library)

Module `RobotProblem` following the 6D template exactly (same factory surface —
`system(; robot_urdf, tstep, domain, Δt_simu, simulator)` even where kwargs are ignored,
`problem()`, `warmup_robot_problem!`, `default_robot_domain()`), so the whole existing
`execution/` infra becomes reusable via a 3-line `elseif "4D"` in
`execution/common/robot_setup.jl`. Contents:

- Geometry + FK/IK ported from Riccardo (harmonize leg lengths with the canonical 6D
  constants: he uses `0.202/0.172` vs `Lthigh = 0.20125`, `Lleg = 0.172` elsewhere — confirm
  with him, default to the canonical values). **Do not trust the ported formulas**: property
  tests are part of the port (see Phase 4).
- **Sound obstacle carving**: OUTER-Lipschitz cell removal, 2-D × 2-D decomposed; ground
  constraint included; returns `UT.set_minus(X, UT.set_union(cells))`.
- **Deduplicated target set** (Set of cell indices).
- **Commensurable-grid constructor**: takes `(hx, τ)`, builds `du = hx/τ`, and *checks* the
  exact-lattice property explicitly; docstring states the bisimulation claim.
- `CENTER_SIMULATION` kept as the approx mode, made legitimate by the lattice check (§3.1) —
  exact, deterministic, one successor per (cell, input), and the cheapest mode available.
- Visualization in the house style: `system_plot!` closure + robot animation
  (port `postproc.jl`), compatible with `animate_trajectory_dashboard`.

### Phase 2 — sampling-time bound (`src/`)

Small documented helper implementing the no-jump theorem (§3.1) — e.g.
`intersample_safe_time_step(hx, U)` — plus a `@warn` in the UGA solver when
`time_step · u_max > hx` while the domain has holes (`SetMinus`). Generic value beyond the
robot; this is the clean, corrected answer to Riccardo's suggestion.

### Phase 3 — input-rate-constrained Dijkstra (`src/optim/discrete_systems/`)

**Is this standard? Yes — it is a classical construction**, which is reassuring:

- In graph terms, constraining consecutive edge labels (`d(u⁻,u) ≤ Δ`) is exactly the
  **turn-restriction / turn-penalty shortest path** problem from road routing, solved since
  Caldwell (1961, *On finding minimum routes in a network with turn penalties*) by running
  Dijkstra on the **line graph** (nodes = edges `(q, u)`), a.k.a. the edge-based or dual
  graph (Winter 2002). Riccardo's per-*node* value table is precisely the known pitfall that
  this literature warns about — labels must live on edges, not nodes.
- In control terms, a bound on `Δu` is the standard **input slew-rate constraint**, handled in
  MPC textbooks (Maciejowski; Rawlings–Mayne) by **state augmentation** `x⁺ = (x, u_prev)` —
  the Δu formulation. Our product automaton is the symbolic-control instance of the same
  augmentation, and the resulting memory-one controller is the textbook outcome.

So we are not inventing an algorithm; we are instantiating a standard one on our automata,
which also extends beyond Riccardo's version (his requires determinism; worst-case Dijkstra on
the product handles non-deterministic automata unchanged).

Correct implementation via the **product automaton** `(q, u_prev)`:

- reuse the existing Dijkstra/value machinery unchanged on the lifted automaton;
- result = a **dynamic controller** (memory = last input), fitting the existing
  `AbstractController` protocol;
- `OptimalControlProblem` stays untouched — the constraint `(d, Δ, target_u)` becomes a
  discrete-solver attribute (e.g. `input_rate_bound`);
- works for non-deterministic automata too (worst-case over successors), unlike his version;
- his 3×3 toy example becomes a unit test.

### Phase 4 — laptop driver + tests

- Driver at `examples/BipedRobot/biped_4d_velocity.jl` (root env — the pure-kinematic model
  needs neither RigidBodyDynamics nor MeshCat): abstraction, obstacle, synthesis, closed-loop
  trajectory, robot animation. Optionally an input-rate-constrained variant once Phase 3 lands.
- Fast tests wired into `test/runtests.jl`: obstacle-carving soundness (a cell grazing the
  preimage must be removed), lattice-exactness check, τ-bound helper, product-Dijkstra on the
  toy automaton, mini end-to-end (small 4-D grid).

**Error vigilance — the port is test-first.** Riccardo's code was never reviewed or tested;
beyond the bugs already found (§3.2), assume more are hiding. Every ported formula gets a
property test *before* the driver relies on it:

- FK/IK round-trip: `get_angular_coordinates ∘ get_cartesian_coordinates = id` on random
  configurations, both stance conventions;
- target set validated *semantically*: for every cell in `compute_target_set`, apply the FK to
  the cell center and check the swing foot lands within tolerance of the requested foothold
  (this cross-checks his hand-derived two-link IK, including the `atan`/matrix-inverse branch);
- obstacle carving validated adversarially: sample points inside carved cells, map through FK,
  assert none hits the (uninflated) obstacle or the ground;
- lattice exactness asserted, then one closed-loop trajectory checked step-by-step against the
  abstract path (bisimulation witnessed numerically);
- code style per the house rules: typed struct fields (his `geometry` struct is untyped and
  lowercase → `Base.@kwdef struct Geometry{T}` with `SVector`-friendly fields), `snake_case`
  functions, docstrings on every exported symbol (docs build errors otherwise), JuliaFormatter
  before each commit.

**Extensibility guardrails**: the obstacle carving is written against a generic
`foot_map(x) -> SVector{2}` + Lipschitz bound rather than hardcoding the biped FK, so the same
helper serves any task-space constraint (and later the moving-obstacle clock product, §5); the
commensurable-grid constructor and the τ-bound helper are model-agnostic from day one.

### Ordering & effort

Phases are independent except 4-depends-on-1 (and the driver variant on 3). Suggested order:
1 → 2 → 4 → 3 (the model + driver give visible results first; the Dijkstra is the largest
new-code item). One commit per phase, `--fast` gate at each phase end.

---

## 5. Beyond the integration (later, with Riccardo)

- **Close the certification gap**: the low-level velocity tracking error `‖v_real − u‖ ≤ ε`
  can be re-imported as a bounded disturbance → `jacobian_bound = 0` stays, growth radius
  becomes `ε·τ`; the abstraction degrades gracefully from bisimulation to alternating
  simulation with a quantified margin. This would make the reduced approach *certified
  end-to-end*.
- Multi-step walking: chain footstep targets (`compute_target_set` at successive footholds)
  with the swap of stance/swing legs — the 6D two-step controller pattern already shows how.
- Moving obstacle: the current carving is static; a time-varying obstacle means product with a
  clock (the UGA `clock`/periodic machinery exists) — natural follow-up once the static case
  is merged.
