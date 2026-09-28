# Underactuated four-link biped: setting up a footstep example in Dionysos

*Session report, 2026-09-28, branch `rb-underactuated-robot`. Nothing below is committed yet.*

## TL;DR

- **Goal.** Build an example in the style of `examples/BipedRobot/biped_4d_velocity.jl`, covering only its `step` scenario. The dynamics come from the new `BipedFourLink` module: the velocity-controlled version, 5-D state, no input accelerations.
- **What exists now.** A problem module, [`problems/BipedRobot/4D_model_underactuated/robot_problem.jl`](#1-the-problem-module), and a runnable driver, [`examples/BipedRobot/biped_4d_underactuated.jl`](#2-the-example-driver). The whole pipeline runs end to end: dynamics → carving → abstraction → controller synthesis → concrete closed-loop simulation.
- **What does not work yet.** No footstep controller was obtained on the real (concrete) system:
  1. **With the obstacle**, the start state is removed from the domain. The sound obstacle carving needs a 7.5 cm safety margin at the affordable resolution, and the swing foot starts closer than that. (The old fully actuated model had the same resolution problem.)
  2. **Without the obstacle**, the abstraction does produce a controller. But on the real dynamics the robot leaves the controller's domain after 7–9 steps. The ankle is an unstable inverted pendulum, and the abstraction rounds the ankle state to cell centres, so it is **not a certified over-approximation** here. The concrete ankle drifts away from what the abstraction predicted. Making the ankle grid 2× finer did not fix it.
- **Cost of trying.** 4 abstraction runs, from 18 k to 5.3 M states. The largest took 14 min and 20 GB of memory.

---

## 0. Background for the reader

Dionysos computes controllers by **abstraction**. The continuous state space is covered by a grid of cells, and each cell becomes a state of a finite automaton. For every (cell, input) pair Dionysos works out which cells the system can reach after one time step `tstep`. A graph algorithm (Dijkstra, for a reach specification) then finds, for every cell, an input that leads to the target. The result is a **state-feedback controller on cells**. It is *concretized* by looking up the cell that contains the real state.

Three Dionysos notions matter for this report:

| Notion | Meaning here |
| --- | --- |
| **`UniformGridAbstraction` solver** | Builds the automaton on a uniform grid. Its `approx_mode` decides how a cell's successors are computed: `GROWTH`/`LINEARIZED` over-approximate the reachable set (sound, certified); `CENTER_SIMULATION` simulates only the cell centre (one successor per input, **not sound in general**). |
| **Carving** (`MP.CellUnion`, `MP.image_blocked_cells`) | Removes cells from the state domain. A cell is removed when some point in it *might* violate a constraint expressed in another space, here "the swing foot touches the obstacle" in Cartesian space. This is made sound by a Lipschitz bound, and the resulting safety margin grows with the cell size. |
| **Concrete closed loop** (`ST.get_closed_loop_trajectory`) | Simulates the *real* dynamics under the synthesized controller. It stops early when the controller has no action for the current cell. When the abstraction is not certified, this is the only real test. |

### The reference example: `biped_4d_velocity.jl` (fully actuated)

The old model is a planar biped with 4 joint angles (stance hip/knee, swing hip/knee), all **actuated in velocity**: the state is the 4 angles, the input is the 4 joint velocities, and the dynamics are `ẋ = u`. That allows a special trick:

- Choosing `tstep · du = dx` makes every input move the state by an exact number of cells. The abstraction is then exact (a bisimulation), so `CENTER_SIMULATION` is *sound for that model* (`MP.is_lattice_exact`).
- The obstacle, an 8 cm × 3 cm triangular bump at the stance foot, is carved out in joint space. `step` asks the swing foot to go from behind the stance foot to a foothold at `(0.2, 0)` m without touching the bump.
- That example already notes that at `dx = 0.1` the carving margin (~5.5 cm) disconnects the free space, so it runs at `dx = 0.05`.

---

## The new robot: `BipedFourLink`

`problems/BipedRobot/4D_model_underactuated/BipedFourLink.jl` (+ `four_link_biped.urdf`) is a **RigidBodyDynamics** model of a planar biped without a torso:

```
            hip  (2 × 0.4 kg motors)
           /   \
   femur  /     \  swing femur        segment lengths:
  0.202  /       \                      tibia 0.172 m, femur 0.202 m
        o knee    o swing knee          (0.4 kg motor at each knee;
        |          \                     links weigh ~20 g, so the robot is
  tibia |           \ swing tibia        essentially point masses at hip/knees)
        |            \
   ankle o (fixed to the ground)   o swing foot
```

| Coordinate | Meaning | Actuated? |
| --- | --- | --- |
| `q₁` | grounded tibia vs. world vertical (the **ankle**) | **no**, it is free |
| `q₂` | grounded femur vs. grounded tibia (stance knee) | yes |
| `q₃` | swing femur vs. grounded femur (hip) | yes |
| `q₄` | swing tibia vs. swing femur (swing knee) | yes |

The key difference from the old model: **the ankle has no motor.** Only `q₂, q₃, q₄` track commanded velocities. The module provides a *velocity-controlled design model* with

```
x = [q̇₁, q₁, q₂, q₃, q₄]     u = [v₂, v₃, v₄]
q̈₁ = −c₁(q, [q̇₁, u]) / M₁₁(q),   q̇₁ = x₁,   (q̇₂, q̇₃, q̇₄) = u
```

This is the first (unactuated) Euler–Lagrange row, with the velocity commands held constant over a step (`input_accelerations = 0`, as you asked). Physically, **the stance leg is an inverted pendulum pivoting on the ankle.** Gravity makes it fall, and the controller can only fight that indirectly, by moving the motor joints and so shifting the centre of mass.

> **State ordering.** Your message wrote the state as `[q₁, q̇₁, q₂, q₃, q₄]`; the module (docstring and `.md`) uses `[q̇₁, q₁, q₂, q₃, q₄]`. I followed the **module's order** so its functions can be called without permuting. Switching would be a small change (the state indices are gathered in `IDX_Q̇1` / `IDX_Q` / `joint_angles`).

---

## What was built

### 1. The problem module

`problems/BipedRobot/4D_model_underactuated/robot_problem.jl`, module `UnderactuatedRobotProblem`, imported as `RP` in the example. Same structure as `4D_model/robot_problem.jl`; it `include`s `BipedFourLink.jl` and uses it as the source of truth:

| Piece | What it does |
| --- | --- |
| `build_model()` | `BipedFourLink.build_model` on the bundled URDF. |
| `Geometry`, `cartesian_coordinates(q, g)` | Forward kinematics of the new URDF (absolute segment angles `q₁`, `q₁+q₂`, `q₁+q₂+q₃`, `q₁+q₂+q₃+q₄`). |
| `swing_foot_gradient_bound` | Per-coordinate Lipschitz bound of the swing-foot position: each joint rotates the rest of the chain, so the bound is the pivot-to-foot distance (ankle 0.748 m, knee 0.576 m, hip 0.374 m, swing knee 0.172 m, `q̇₁` 0). |
| **`VelocityControlledDynamics(model)`** | Callable `(x, u) -> ẋ`. Same formula as `BipedFourLink.velocity_controlled_derivative`, but it reuses **preallocated per-thread** RigidBodyDynamics buffers on the module's own mechanism and returns an `SVector`. See [step 3](#step-3--speed-of-the-modules-dynamics-too-slow--fast-path). |
| **`check_against_module(model)`** | Compares the fast dynamics and the forward kinematics with `BipedFourLink` at 50 sample points and errors on any mismatch. The example calls it at startup. |
| `infeasible_cells` | Cells to remove from the domain: swing foot possibly in the obstacle (sound, `MP.image_blocked_cells`), or surely below the ground. It is computed on the 4-D angle grid, because the foot doesn't depend on `q̇₁`, and then lifted along the `q̇₁` axis. |
| `target_cells` / `target_set` | Cells where the swing foot is on the foothold. `(q₁, q₂)` are swept over the grid; `(q₃, q₄)` come from two-link inverse kinematics with the knee bent backwards (`q₄ ≥ 0`); any `q̇₁` is accepted. |
| `balanced_start(f; φ3, q2, q4)` | A rest posture whose ankle is in (unstable) static equilibrium, `q̈₁ = 0`, found by bisection on `q₁`. See [step 4](#step-4--open-loop-behaviour-the-old-start-pose-is-unusable). |
| `default_discretization(; dx, dq1, dq̇1, tstep)` | 5-D state grid `(dq̇1, dq1, dx, dx, dx)`, 3-D input grid with `du = dx / tstep`, so the three motor coordinates stay on the lattice. |
| `system(; model, domain, removed_cells)` | A `ConstrainedBlackBoxControlContinuousSystem`. The solver discretizes it with RK4 at its `time_step`. |
| `step_problem`, `draw_robot!`, `system_plot!` | Reach specification and plotting, same interface as the old module. |

### 2. The example driver

`examples/BipedRobot/biped_4d_underactuated.jl` follows the old example's layout: scenario → carving → abstraction → footstep synthesis → concrete closed loop → plots. It differs as follows:

- **Only the `step` scenario**: the same bump `_triangle(0.0, 0.04, 0.03)` and foothold `(0.2, 0)`.
- **No slew-rate pass**. It belongs with the `input_accelerations` version you want to add later.
- **`approx_mode = CENTER_SIMULATION`**, the same mode the existing 6D/8D RigidBodyDynamics biped models use. A growth bound would need a Jacobian bound of `q̈₁`, which isn't available for a RigidBodyDynamics black box. The file header says clearly that the result is **not certified**.
- Prints the carving margin, and **errors with an explanation if the start cell was carved out**.
- Prints why the concrete run stopped ("outside the controller's domain") and the range the ankle covered.
- Extra plot: the ankle `q₁, q̇₁` against time.
- Same override-by-global pattern as the old example: `DISC`, `OBSTACLE`, `X_BOUNDS`, `ONE_JOINT_PER_STEP` (default `false`, i.e. all 27 input combinations), `RUN_PLOTS`.

Run command, as written in the file header:

```bash
julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_underactuated.jl
```

### 3. Environment

The other biped examples run with `--project=test`, but that environment has neither `RigidBodyDynamics` nor `MeshCatMechanisms`, and `BipedFourLink` needs both. `problems/BipedRobot/Project.toml` has both but lacked `LazySets` and `Plots`, so **I added `LazySets` and `Plots` to it**. Its `Manifest.toml` is gitignored. All experiments ran in a scratch environment with the same packages.

---

## Step-by-step log (including what failed)

### Step 1: understanding the existing code

I read the old example, the old problem module, `BipedFourLink.jl` / `.md`, the URDF, and how the 6D/8D RigidBodyDynamics bipeds are abstracted (they use `CENTER_SIMULATION`). I also read the `UniformGridAbstraction` approximation modes and the carving/lattice helpers. The helpers `MP.image_blocked_cells` and `MP.cells_where` work in any dimension and can be reused. `MP.swept_input_filter`, the sound inter-sample check, assumes `ẋ = u`, so it **does not apply** to the ankle.

### Step 2: verifying the kinematics ✅

I derived the forward kinematics from the URDF and compared them with the module's `foot_position` / `world_position` at random postures. They agreed to machine precision. Gravity sign check: `q₁ = 0.1` at rest gives `q̈₁ = +2.3 rad/s²`, so the ankle falls further, as an inverted pendulum should (`ω² ≈ 23 s⁻²`, time constant ≈ 0.2 s).

### Step 3: speed of the module's dynamics (too slow) → fast path ✅

`velocity_controlled_derivative` takes **174 µs and 7,180 allocations per call**. It builds a new mechanism state and dynamics result every time and computes the mass matrix and the bias forces in two separate passes. One RK4 step with 5 substeps costs 20 calls, about 3.5 ms. A few million (cell, input) pairs would take hours.

`VelocityControlledDynamics` computes the same quantities, `M₁₁` via `mass_matrix!` and `c₁` via `dynamics_bias!`, into per-thread buffers:

| | time / call | allocations |
| --- | --- | --- |
| `BipedFourLink.velocity_controlled_derivative` | 174 µs | 7,180 |
| `VelocityControlledDynamics` | **1.6 µs** | **0** |
| one RK4 step (5 substeps) with the fast path | 33 µs | 0 |

It agrees with the module at 200 random `(x, u)` points to 1e-9. That check ships as `check_against_module`.

### Step 4: open-loop behaviour, the old start pose is unusable ⚠️

The old example starts with the hip behind the stance foot. In new coordinates that is roughly `q = (−0.2, 0, 0.4, 0)`. With the joints frozen, the ankle there goes from `q₁ = −0.2` to `−0.53` in 0.3 s and `−5.2` in 1 s: the robot falls backwards at once. So I added `balanced_start`. It keeps the swing leg straight and 0.2 rad behind the hip, as in the old start, but picks `q₁` so the ankle is in static equilibrium. The result is snapped to the nearest cell centre, so the concrete run starts exactly where the abstraction does.

### Step 5: run 1, coarse smoke test ❌ (no controller)

`dx = 0.2, dq1 = 0.1, dq̇1 = 0.2`, first bounds `|q̇₁| ≤ 1, |q₁| ≤ 0.3, |q₂| ≤ 0.8, |q₃|, |q₄| ≤ 1.2`. That gave 18,090 states and 70,852 transitions, built in ~3 s. No controller, 0 steps. The pipeline itself worked. In hindsight the start cell was almost certainly carved out here as well: at this resolution the margin is ~15 cm. See step 8.

### Step 6: is the step physically possible within the bounds? ⚠️

For every posture that puts the swing foot on the foothold, I computed `q̈₁` at rest:

- With `|q₁| ≤ 0.3, |q₂| ≤ 0.8`, **every target posture falls forward at `q̈₁ ≥ 4.5 rad/s²`**.
- The only balanced target postures need `q₁ ≈ −0.7, q₂ ≈ +1.2` (a deep, backward-bent stance knee).

So the footstep has to be a **controlled fall into the foothold**, and the ankle velocity range must leave room for it. I widened the bounds to `q̇₁ ∈ [−1.5, 1.5]`, `q₁ ∈ [−0.4, 0.4]`, `q₂ ∈ [−0.8, 0.8]`. I narrowed the swing leg to its useful side: `q₃ ∈ [−1.2, 0.4]` (swing forwards) and `q₄ ∈ [−0.2, 1.2]` (knee bends backwards).

### Step 7: run 2, medium resolution with the obstacle ❌ (no controller)

`dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1`, all 27 inputs. That gave **716,909 states, 19.4 M forward images, 10.9 M transitions, built in 113 s, 3.4 GB peak**. Dijkstra found no controller.

Two failed diagnostic attempts followed, both bugs in my throwaway scripts. First, the implicit mapping's state ids are not `1:n_states`, so an array indexed by id overflowed. Second, the id lookup needs `Int`, not `Int32`. Each cost a full rebuild. After that the scripts cached the transitions to disk before analysing them.

### Step 8: diagnosis, the start cell is carved out ❌

The fixed diagnostic showed that **the start cell has no outgoing transitions at all.** The carving had removed it. The margin is

```
margin = Σᵢ (Lipschitz bound)ᵢ · hᵢ / 2
       = 0.748·dq1/2 + (0.576 + 0.374 + 0.172)·dx/2  =  7.5 cm   (dq1 = 0.05, dx = 0.1)
```

At the start, the swing foot is about 7 cm from the bump, inside that margin. This is the same resolution problem the old example documents, but worse, because here the unactuated ankle swings the *whole robot* (0.748 m lever arm) and the stance knee swings everything above it. I added an explicit check to the example: it now **errors with this explanation** instead of silently finding nothing. **With the default settings the example therefore stops at this check.**

### Step 9: run 3, medium resolution without the obstacle ⚠️ (abstract ✅, concrete ❌)

Same grid, `OBSTACLE = nothing`. With only below-ground cells carved, that gave 1,267,532 states, 34.2 M forward images and 21.7 M transitions, built in 197 s. **Dijkstra succeeded**: 9 abstract steps.

Automaton analysis (the abstraction is deterministic, so "a controller exists" is the same as "a path exists"):
- 292,399 states are reachable from the start (max depth 38).
- 340,757 states can reach the target, 14,207 of them at rest (`q̇₁ = 0`).
- 3,877 target states.

**Concrete closed loop: 9 steps, then the controller has no action for the current cell.** The ankle reached `q₁ ≈ 0.11`, `q̇₁ ≈ 1.1`. The real ankle had drifted into a cell that is not in the controller's domain.

### Step 10: run 4, finer ankle grid without the obstacle ⚠️ (abstract ✅, concrete ❌)

Hypothesis: the drift comes from rounding the ankle state. I halved both ankle steps (`dq1 = 0.025`, `dq̇1 = 0.05`). That gave **5,313,186 states, 143.5 M forward images, 91.0 M transitions, built in 844 s, Dijkstra 31 s, 20.2 GB peak**, with a 6.5 cm margin. **Dijkstra succeeded again.** Concrete closed loop: **7 steps**, stopping outside the controller's domain at `x = [0.94, 0.22, −0.30, −0.40, 0.10]`, with the robot falling forward. The finer ankle grid did **not** fix the mismatch. Following your instruction, I stopped experimenting here.

### Summary of runs

| Run | grid `(dq̇1, dq1, dx)` | obstacle | states | transitions | build | peak mem | abstract | concrete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (0.2, 0.1, 0.2) | yes | 18 k | 71 k | ~3 s | — | ✗ | — (start carved, in hindsight) |
| 2 | (0.1, 0.05, 0.1) | yes | 717 k | 10.9 M | 113 s | 3.4 GB | ✗ | — (start carved) |
| 3 | (0.1, 0.05, 0.1) | no | 1.27 M | 21.7 M | 197 s | — | ✓ 9 steps | ✗ stops after 9 steps |
| 4 | (0.05, 0.025, 0.1) | no | 5.3 M | 91.0 M | 844 s | 20.2 GB | ✓ | ✗ stops after 7 steps |


---

## Why it doesn't work (yet)

1. **The obstacle needs a finer motor grid than is affordable in 5-D.** The margin is dominated by the ankle and stance-knee terms. Getting close to the old model's working margin (2.7 cm) needs about `dx = 0.05` and `dq1 = 0.025`. With the current bounds that is ~32 M cells × 27 inputs ≈ 875 M RK4 roll-outs: roughly 50 min on 10 threads and well over 20 GB.
2. **`CENTER_SIMULATION` on an unstable, nonlinear coordinate isn't robust.** The motor joints are still exact (lattice), but `(q₁, q̇₁)` are rounded to cell centres at every step. The abstraction then predicts one successor per input, and the synthesized controller is only defined along those predicted paths. The real, unrounded ankle drifts off them within a few steps, the more so because it diverges exponentially (time constant ≈ 0.2 s).
3. **The velocity alphabet is slow.** The lattice trick fixes the joint speed at `du = dx / tstep = 1 rad/s`. The controller must keep the ankle balanced while the swing leg travels ~0.5 rad and flexes its knee over the bump. Every target posture is falling forward (step 6), so the last part of the step is a race against the fall.
4. **Inter-sample safety is not checked.** Ankle moves span several cells per step (the solver's warning reports up to 16), and the swept-cell check only applies to `ẋ = u`. This doesn't cause the failures above, but it is another reason the result is not certified.

## Not done / not verified

- **Plotting path untested.** No run produced a trajectory worth drawing, so the dashboard/PNG code has not run. It is guarded to skip when the trajectory is empty.
- **The documented run command** (`--project=problems/BipedRobot`) was not tested end to end. All runs used a scratch environment with the same packages. The environment instantiated without errors afterwards.
- **No unit tests.** `test/Project.toml` lacks RigidBodyDynamics. The fast dynamics are validated at runtime by `check_against_module`, not by the test suite.
- **Formatting.** My two new files are JuliaFormatter-clean. `BipedFourLink.jl` (your file) is not in the repo style, e.g. functions without explicit `return`. CI runs `format(".")` on the whole repo, so the format check will fail until it is formatted.
- **Slew-rate / `input_accelerations` version**: not started.
- **Not committed**: the new `robot_problem.jl`, the example, and the `Project.toml` change.

## Options for next steps

Roughly from cheapest to most involved:

1. **Make the concrete loop more robust without certification.** Use `RANDOM_SIMULATION` with several samples per cell. The automaton becomes non-deterministic, so the controller must work for the spread of real successors, not just the centre's. It costs ×`n_samples` in build time and may make the problem infeasible, but that would be an honest answer.
2. **Tighter carving.** Replace the global Lipschitz box bound with a per-cell bound, for example the actual Jacobian at the cell centre plus a second-order term. The margin (7.5 cm) is what blocks the obstacle scenario at `dx = 0.1`. This would touch the library's carving helper or add a model-specific variant.
3. **Relax the scenario.** Use a smaller bump, start the swing foot further back, use a shorter foothold (so target postures fall less), or restrict `q₂` to the human-like side to shrink the grid.
4. **Faster joints.** Two speed levels per joint (`speed_levels = 2`, as in the old example) mean 5³ = 125 inputs. That is costly, and the multi-cell moves have no sound swept check here.
5. **Certified route.** A `GROWTH` abstraction needs a bound on the Jacobian of `q̈₁` over each cell. It could come from interval arithmetic on the RigidBodyDynamics terms (sound), or from sampling with ForwardDiff (not sound) as a first step.
6. **The `input_accelerations` model.** It is more faithful at input switches, and it pairs naturally with the slew-rate (`BoundedInputVariation`) pass of the old example.

## Earlier in the session (unrelated)

I added `JuliaLint.toml`, which turns off the Julia VS Code extension's "computed include" warning outside `src/`/`ext/`. It is committed as `703cdb9d3` on this branch.

## Reproducing

```bash
# default: step scenario at dx = 0.1 → stops at the "start cell carved out" check
julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_underactuated.jl
```

Run 3 (no obstacle, medium grid), from a Julia session in that environment:

```julia
OBSTACLE = nothing
include("examples/BipedRobot/biped_4d_underactuated.jl")
```

Run 4 additionally needs, before the `include`:

```julia
include("problems/BipedRobot/4D_model_underactuated/robot_problem.jl")
DISC = UnderactuatedRobotProblem.default_discretization(; dx = 0.1, dq1 = 0.025, dq̇1 = 0.05, tstep = 0.1)
```
