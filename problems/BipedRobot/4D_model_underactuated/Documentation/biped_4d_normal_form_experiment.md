# The `biped_4d_normal_form.jl` experiment

`examples/BipedRobot/biped_4d_normal_form.jl` computes a certified abstraction of an underactuated
four-link biped and asks Dionysos for a controller that performs **one footstep over an obstacle**.
This document explains what it does, lists every parameter you can tune (with its default), and
walks through the code step by step. The modules it relies on are documented in
[`symbolic_multibody_and_local_abstraction.md`](symbolic_multibody_and_local_abstraction.md) (the
file `robot_problem.jl` in its §4).

**Contents**

1. [What the experiment does](#1-what-the-experiment-does)
2. [Running it, and what a default run prints](#2-running-it-and-what-a-default-run-prints)
3. [Tunable parameters](#3-tunable-parameters)
4. [Step-by-step walk-through of the code](#4-step-by-step-walk-through-of-the-code)
5. [Diagnosing a failed run](#5-diagnosing-a-failed-run)

---

## 1. What the experiment does

### 1.1 The robot

A planar biped with four joint angles and **no torque at the ankle**. In the notation of the code
(angles counter-clockwise, `q = 0` is the robot standing straight with both legs side by side):

| Coordinate | Meaning | Actuated? |
| :--- | :--- | :--- |
| `q₁` | grounded ankle, relative to the vertical | **no** |
| `q₂` | grounded knee (stance femur relative to stance tibia) | yes |
| `q₃` | hip (swing femur relative to stance femur) | yes |
| `q₄` | swing knee (swing tibia relative to swing femur) | yes |

The three motor joints are commanded **in velocity**: `u = (v₂, v₃, v₄)` is held constant over each
sampling period. The ankle obeys the first Euler–Lagrange row of the multibody model: it falls like
an inverted pendulum, and the controller can only influence it through the posture of the motor
joints (inertial coupling and gravity torque).

### 1.2 The model: normal form with the ankle *momentum*

The multibody model is built symbolically with `Multibody2D`, brought to normal form
`x = (q, v_a, p_u)` and reduced to a 5-dimensional model in which the unactuated ankle is
described by its **generalized momentum** `p₁ = ∂L/∂q̇₁ = (D q̇)₁` instead of its velocity:

```
state   x = [p₁, q₁, q₂, q₃, q₄]          input  u = [v₂, v₃, v₄]
ẋ = [ ṗ₁, q̇₁, v₂, v₃, v₄ ],     q̇₁ = (p₁ − D₁ₐ v_a) / D₁₁ ,     ṗ₁ = ∂L/∂q₁ .
```

`p₁` takes the slot of the ankle velocity of the original example (`biped_4d_underactuated.jl`), so
the geometry, obstacle carving, target set and plots of `robot_problem.jl` apply unchanged. Because
`p₁ ≈ D₁₁ q̇₁ + D₁ₐ·v_a` with `D₁₁ ≈ 0.14 kg m²`, a momentum of `0.02 N m s` corresponds to about
`0.15 rad/s` of ankle velocity.

### 1.3 The abstraction: local growth bounds

The state space is divided into hyper-rectangular cells and the inputs into a finite set of velocity
commands. For every `(cell, input)` pair the **reachable set after one sampling period** is
over-approximated by a box (the cell centre is simulated, and the box radius comes from a growth
bound computed *locally*, with interval arithmetic, for that cell and input; see the module
document). The result is a nondeterministic finite transition system on which a controller is
synthesized against the worst case, so a controller found on the abstraction is guaranteed to work on
the concrete system (up to numerical integration and the boundary convention of §3.5).

### 1.4 The specification

* **Foot-height constraint.** The swing foot must stay above the ground `z = 0`. Cells where the
  swing foot is *surely* below the ground are removed from the state space.
* **Obstacle at the origin.** A triangular bump (8 cm wide, 3 cm high) sits at the stance foot.
  Every cell whose swing foot *might* touch it is removed, which forces the foot to lift.
* **Footstep.** Starting from a configuration `q0`, with **any ankle momentum** in a given range, reach
  a configuration placing the swing foot on a foothold 20 cm in front of the stance foot, again
  with any ankle momentum in a given range. The default ranges are the whole interval of definition
  of the state space, because the swing foot has no ground reaction force and requiring rest at the
  start or at the end would be needlessly conservative.

Both carvings are *sound*: they use a Lipschitz bound of the swing-foot map over a cell (the
"Lipschitz carving margin", 7.5 cm with the default grid), so they remove more than strictly
necessary.

---

## 2. Running it, and what a default run prints

### 2.1 Requirements

An environment with Dionysos (developed from this repository: the example locates
`problems/BipedRobot/4D_model_underactuated/` through `pathof(Dionysos)`), Symbolics,
IntervalArithmetic, MathematicalSystems, LazySets, StaticArrays, MathOptInterface and Plots (loaded by
`robot_problem.jl`). RigidBodyDynamics is **not** needed. The `problems/BipedRobot` project provides
most of them; Symbolics and IntervalArithmetic must be added.

```
julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_normal_form.jl
```

Use several threads (`-t auto`): the abstraction is the expensive step and is threaded.

### 2.2 A default run (10 threads)

```
Symbolic normal form: NormalForm(8 states, 3 inputs; actuated = [2, 3, 4], unactuated = [1])
Carving the obstacle and the ground out of the state domain…
Lipschitz carving margin: 7.5 cm
x0 = [0.0, -0.1, 0.0, -0.5, 0.0]  (ṗ₁ at x0: -0.0412)
swing foot at x0: [-0.1738, 0.0635]
[ Info: Number of states: 176225
[ Info: Number of inputs: 27
[ Info: Number of forward images: 4758075
transitions: 18820992
local growth bound statistics: (calls = 4758075, enclosure_failures = 0, unbounded = 0, improved = 4758075)
Optimal control terminated with success = false
footstep success: false
steps: 0, reached: false
```

| Step | Time |
| :--- | :--- |
| symbolic normal form (mostly compilation) | ≈ 11 s |
| carving the obstacle and the ground | ≈ 4.6 s |
| `LocalGrowthBound` (symbolic Jacobian, compilation) | ≈ 5 s |
| abstraction (4.76 M `(cell, input)` pairs) | ≈ 112 s |
| controller synthesis | ≈ 1.4 s |

**With the defaults the synthesis finds no controller** (`footstep success: false`): the
closed-loop run has no action at `x0` (`steps: 0`) and the plots are skipped. The pipeline itself
works end to end: with a narrower initial range or a finer grid the problem is easier (§5). A
transitions warning ("Some transitions jump 7.0 cells per step while the state domain has holes")
is printed by Dionysos because the obstacle leaves holes in the domain and the abstraction does not
sweep the trajectories between samples.

---

## 3. Tunable parameters

Parameters come in four kinds, from the easiest to change to the most invasive.

### 3.1 Variables you define *before* including the file

The example reads these with `@isdefined`, so you can run it from a driver script without editing
it. **Defaults are used for any variable you do not define.**

```julia
# driver.jl
using StaticArrays, Dionysos
const MP = Dionysos.Mapping
dx, dq1, dp, tstep = 0.1, 0.1, 0.06, 0.05                  # a coarser grid
du = dx / tstep
DISC = (; state_grid = MP.GridFree(zeros(SVector{5, Float64}), SVector(dp, dq1, dx, dx, dx)),
          input_grid = MP.GridFree(zeros(SVector{3, Float64}), SVector(du, du, du)),
          tstep, u_max = du, du)
INITIAL_VELOCITY = (-0.05, 0.05)                            # start with a small ankle momentum
RUN_PLOTS = false
include(joinpath(dirname(dirname(pathof(Dionysos))), "examples", "BipedRobot", "biped_4d_normal_form.jl"))
```

(`include` resolves a relative path against the directory of the *script that calls it*, so use an
absolute path such as the one above, derived from where Dionysos is installed.)

| Variable | Default | Meaning and constraints |
| :--- | :--- | :--- |
| **`DISC`** | `RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05)` | Grids and sampling time, a `NamedTuple` with fields `state_grid`, `input_grid`, `tstep`, `u_max` (see §3.2). |
| **`X_BOUNDS`** | `(x_lb, x_ub) = (SVector(-0.2, -0.3, -0.5, -0.8, -1.0), SVector(0.2, 0.2, 0.5, 1.0, 0.2))` | State box `[p₁, q₁, q₂, q₃, q₄]` (§3.3). |
| **`Q0`** | `SVector(-0.1, 0.0, -0.5, 0.0)` | Initial **configuration** `[q₁, q₂, q₃, q₄]` (snapped to the grid). Must be in the domain and not carved out (§3.4). |
| **`INITIAL_VELOCITY`** | `nothing` = the whole interval `[x_lb[1], x_ub[1]] = [-0.2, 0.2]` | `(vlo, vhi)` range of the first state coordinate `p₁` at the start. The controller must work from every velocity cell the range touches. |
| **`TARGET_VELOCITY`** | `nothing` = the whole interval `[-0.2, 0.2]` | `(vlo, vhi)` range of `p₁` at the end. The target is made of the cells lying entirely inside the range. |
| **`OBSTACLE`** | `[_triangle(0.0, 0.04, 0.03)]` | A vector of `LazySets` polygons (a triangle of half-width 4 cm and height 3 cm centred at `x = 0`), or `nothing` for no obstacle (only the ground constraint). |
| **`ONE_JOINT_PER_STEP`** | `false` | `true` allows at most one moving motor joint per step (7 inputs instead of 27, a much smaller automaton). |
| **`RUN_PLOTS`** | `true` | `false` skips the figures and the animation. |

### 3.2 Grid cell sizes, sampling time and input set: `DISC`

`RP.default_discretization(; dx, dq1, dq̇1, tstep)` builds

* the **state grid** `GridFree(origin = 0, steps = (dq̇1, dq1, dx, dx, dx))`: steps of
  `(p₁, q₁, q₂, q₃, q₄)`. In this experiment `dq̇1` is the step of the **momentum** `p₁`;
* the **input grid** `GridFree(origin = 0, steps = (du, du, du))` with `du = dx / tstep`;
* the sampling time `tstep` and `u_max = du`.

| Parameter | Default | Meaning | Effect and constraints |
| :--- | :--- | :--- | :--- |
| `dx` | `0.1` rad | cell size of `q₂, q₃, q₄` | with `tstep·du = dx` every command moves a motor joint by *exactly one cell* per step (the lattice trick): the three motor coordinates are then exact and the abstraction only has uncertainty in `(q₁, p₁)` |
| `dq1` | `0.05` rad | cell size of the ankle angle `q₁` | finer = tighter ankle model, more cells; also sets the carving margin |
| `dq̇1` | `0.02` N m s | cell size of the ankle momentum `p₁` | finer = better, more cells (`0.02 ≈ 0.15 rad/s`) |
| `tstep` | `0.05` s | sampling period (input held constant) | must be short enough for the flow enclosure to exist: `tstep · ‖∂f/∂x‖ ≲ 1`, and here `‖∂f/∂x‖ ≈ 9`; `0.1 s` is marginal. Shorter steps also make `du` larger |
| `du = u_max` | `dx/tstep = 2.0` rad/s | commanded joint speed | the input set is `{-du, 0, du}³` (27 combinations); `U = [-u_max, u_max]³` |

A custom `DISC` needs `state_grid`, `input_grid`, `tstep` and `u_max`. Keep `tstep · du = dx` if you want the exact lattice.

### 3.3 State-space boundaries: `X_BOUNDS`

`X_BOUNDS = (x_lb, x_ub)`, each an `SVector{5}` in the order `[p₁, q₁, q₂, q₃, q₄]`:

| Coordinate | Default `x_lb` | Default `x_ub` | Why |
| :--- | ---: | ---: | :--- |
| `p₁` ankle momentum | `-0.2` | `0.2` | covers ankle velocities up to about 1.2 rad/s in magnitude (`p₁ ≈ 0.14 q̇₁ + D₁ₐ v_a`) |
| `q₁` ankle angle | `-0.3` | `0.2` | the robot falls *forwards*, towards negative `q₁`, so there is more room on that side |
| `q₂` grounded knee | `-0.5` | `0.5` | |
| `q₃` hip | `-0.8` | `1.0` | the swing hip swings forwards for `q₃ > 0`; the start is at `q₃ = -0.5` |
| `q₄` swing knee | `-1.0` | `0.2` | the knee bends backwards for `q₄ < 0` |

Only cells lying **entirely** inside the box are kept (Dionysos' `INNER` mode), so a box whose faces
are not on a cell boundary loses a partial cell: with the default grid the box contains
**19 × 9 × 9 × 17 × 11 = 287 793 cells** (not the naive `width / step` = 432 000), and **176 225**
remain after carving the obstacle and the ground. A wider box such as
`(SVector(-0.3, -0.4, -0.8, -0.8, -1.2), SVector(0.3, 0.4, 0.8, 1.2, 0.2))` has about three times more
cells and takes proportionally longer.

### 3.4 Initial and target sets

| Item | Default | Defined where | Meaning |
| :--- | :--- | :--- | :--- |
| initial configuration | `Q0 = [q₁,q₂,q₃,q₄] = [-0.1, 0, -0.5, 0]` | `Q0` | stance leg leaning 0.1 rad forwards, knee straight, the **swing leg straight and 0.6 rad behind the hip** (`q₁ + q₃`), which puts the swing foot at `x ≈ -0.17 m`, 13 cm behind the obstacle's edge |
| initial velocity (`p₁`) | whole interval `[-0.2, 0.2]` | `INITIAL_VELOCITY` | any ankle momentum |
| initial set shape | box of half-width `0.25 h` around the centre of the cell of `Q0`, times the velocity range | `RP.step_problem` (fixed) | the controller must work from every velocity cell the range touches |
| foothold | `SVector(0.2, 0.0)` m | `foothold` in the file | the swing foot must end 20 cm in front of the stance foot, on the ground |
| target velocity (`p₁`) | whole interval `[-0.2, 0.2]` | `TARGET_VELOCITY` | any ankle momentum at the end |
| target configurations | all `(q₁, q₂)` cells, `(q₃, q₄)` by inverse kinematics with the knee bent backwards | `RP.target_cells(…; nsub = 3)` (fixed) | cells containing a configuration with the stance foot at the origin and the swing foot at the foothold; `nsub = 3` samples per axis per cell |

Constraints to keep in mind:

* the cell of `Q0` must be in the domain and **not carved out**: the swing foot at `Q0` must be farther
  than the **Lipschitz carving margin** from the obstacle (7.5 cm with the default grid; it scales
  with the cell sizes) and above the ground. The script stops with an error otherwise;
* the foothold must be reachable with `(q₃, q₄)` inside `X_BOUNDS`, otherwise `target_set` reports "no
  reachable target configuration";
* velocity ranges outside the domain interval are rejected, reversed ranges too, and a target
  range must contain at least one whole cell (`≥ dq̇1` wide).

### 3.5 Quantities fixed in the code (edit the file to change them)

These have no override variable; each is one line in `biped_4d_normal_form.jl` unless said otherwise.

| Item | Default | Where | Meaning |
| :--- | :--- | :--- | :--- |
| physical parameters | `L1 = 0.202`, `L2 = 0.172` m, `m1 = 0.0202`, `m2 = 0.0172` kg, `mM = 0.4` kg, `g = 9.81` | `values = RP.default_parameters()` | femur / tibia length, femur / tibia mass, motor mass, gravity. The URDF values; the kinematics of the carving (`RP.Geometry(model, values)`) is compiled from the same symbolic model and `values`, so dynamics and constraints cannot disagree |
| inertia model | uniform thin-rod links `m L²/12`; motors point masses (`I = 0`); hip `2 mM`; swing foot massless | `RP.symbolic_model()` | the model structure (frames, bodies): the robot of `four_link_biped.urdf` |
| scenario name | `:step_normal_form` | `scenario_name` | only used in the output file names |
| `LocalGrowthBound` options | `refinements = 1`, `ngrowthbound = 5`, `nsystem = 5`, `max_enclosure_iterations = 10`, `enclosure_inflation = 0.1`, `shrink_enclosure = false`, `boundary_tolerance = 1e-9` | the call `LGA.LocalGrowthBound(f_sym, x_sym, u_sym, X, U_box, disc.tstep)` | add them as keywords; see below |
| abstraction options | `print_level = 1`, `use_implicit_mapping = true`, `execution_backend = SY.ThreadedBackend()` | the call `LGA.build_abstraction(…)` | verbosity, implicit state mapping (low memory), threading; any other optimizer attribute can be added as a keyword |
| closed-loop horizon | `400` steps | `ST.get_closed_loop_trajectory(…, x0, 400; …)` | maximum number of steps of the concrete run |
| closed-loop start | `q0` with `p₁` at the middle of the initial range, snapped to the grid | `v0`, `x0` | the concrete run starts here (it does not sample the initial set) |
| plot window and files | `xlims = (-0.6, 0.6)`, `ylims = (-0.05, 0.55)`, `snapshot_every = 5` (defaults of `RP.plot_footstep`); `…_step_normal_form.png / .gif / _ankle.png` next to the script | the call `RP.plot_footstep(…)` in the visualization section | add them as keywords |
| ground constraint | foot below `z = -margin` is removed | `RP.infeasible_cells` | the margin makes the ground constraint tolerant so that footholds on the ground stay |

The `LocalGrowthBound` options that matter:

| Option | Effect |
| :--- | :--- |
| `refinements` | `0` = the plain global growth bound of Dionysos (use it as a baseline); `1` = one local refinement (default); more can tighten further |
| `max_enclosure_iterations`, `enclosure_inflation` | how hard the flow enclosure is searched before a cell falls back to the global bound |
| `shrink_enclosure` | `true` uses a tighter verified box for the Jacobian bound (slightly tighter images) |
| `boundary_tolerance` | fraction of the cell radius by which images are deflated so that an image that only *touches* a neighbour does not overlap it. With the exact lattice of the motor joints, `0` gives 3 successors per motor joint instead of 1 and synthesis becomes hopeless |

### 3.6 Cost of a change

* **Number of `(cell, input)` pairs** = (cells after carving) × (number of inputs = 27, or 7 with
  `ONE_JOINT_PER_STEP`). Default: 176 225 × 27 = 4 758 075.
* **Time of the abstraction** ≈ 24 µs per pair of wall time on 10 threads (≈ 0.24 ms of CPU per pair
  of interval arithmetic), i.e. about 2 minutes by default. It is linear in the number of pairs.
* **Memory** grows with the number of transitions (18.8 M by default). Halving every ankle step
  multiplies cells and transitions by about 4.
* Halving `dx` multiplies the cells by 8 (three motor axes). It also halves the commanded joint speed
  `du = dx / tstep` unless `tstep` is halved too (which doubles the number of steps of a footstep).

---

## 4. Step-by-step walk-through of the code

The file is organised in five numbered sections; the paragraphs below follow it.

### 4.0 Header, packages and includes

```julia
using StaticArrays, MathematicalSystems, Dionysos
import MathOptInterface as MOI, LazySets, Plots, Symbolics
const DI = Dionysos; const UT = DI.Utils; const ST = DI.System; const PR = DI.Problem
const MP = DI.Mapping; const SY = DI.Symbolic; const OP = DI.Optim; const AB = OP.Abstraction
```

The standard Dionysos aliases (`UT` utilities, `ST` systems, `PR` problems, `MP` mappings/grids,
`SY` symbolic models, …).

```julia
const PROBLEM_DIR = joinpath(dirname(dirname(pathof(Dionysos))), "problems", "BipedRobot", "4D_model_underactuated")
include(joinpath(PROBLEM_DIR, "robot_problem.jl"));          import .UnderactuatedRobotProblem as RP
include(joinpath(PROBLEM_DIR, "LocalGrowthAbstraction.jl")); import .LocalGrowthAbstraction as LGA
```

* `RP` (`robot_problem.jl`) provides the **symbolic model** (`symbolic_model`, `reduced_dynamics`,
  `substitute_parameters`; it includes `Multibody2D.jl`), the **geometry** (`Geometry`,
  `cartesian_coordinates`), the **constraints** (`infeasible_cells`), the **target set**
  (`target_set`), the **problem** (`step_problem`) and the **plots**. It does not depend on
  RigidBodyDynamics.
* `LGA` (`LocalGrowthAbstraction.jl`) provides the local growth bound and the abstraction builder.

### 4.1 Section 1: symbolic model → numerical dynamics

```julia
values  = RP.default_parameters()
@time model = RP.symbolic_model()
reduced = RP.reduced_dynamics(model)
f_sym   = RP.substitute_parameters(model, reduced.f, values)
x_sym, u_sym = reduced.x, reduced.u
geometry = RP.Geometry(model, values)
```

1. `values` are the numbers (§3.5).
2. `symbolic_model()` builds the four frames and eight points of the biped with **symbolic**
   parameters and computes its normal form (8 states `(q₁…q₄, v₂, v₃, v₄, p₁)`, 3 inputs). About
   11 s, almost all of it compilation.
3. `reduced_dynamics(model)` keeps the `p₁` and `q̇` rows: a 5-state, 3-input velocity-controlled
   model with the ankle momentum **first**, `x = (p₁, q₁…q₄)` and `ẋ = (ṗ₁, q̇₁, …, q̇₄)`, in the slot
   of the ankle velocity of the original example, `u = (v₂, v₃, v₄)`. Both `Multibody2D` and
   `robot_problem.jl` use the counter-clockwise convention, so nothing else is adapted. The
   parameters are still symbolic.
4. `substitute_parameters` puts `values` in. What remains is symbolic in `x` and `u` only.
5. `RP.Geometry(model, values)` holds the two segment lengths and the forward kinematics of the
   joints, **compiled from the symbolic positions of the same `model`** (grounded knee, hip, swing
   knee, swing foot). The carving and the target set computed below with `RP` are therefore valid for
   the symbolic model by construction: there is no second copy of the kinematics to check.

### 4.2 Section 2: scenario, grids, domain, start and carving

```julia
obstacle = @isdefined(OBSTACLE) ? OBSTACLE : [_triangle(0.0, 0.04, 0.03)]
foothold = SVector(0.2, 0.0)
disc = @isdefined(DISC) ? DISC : RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05)
x_lb, x_ub = @isdefined(X_BOUNDS) ? X_BOUNDS : (SVector(-0.2, -0.3, -0.5, -0.8, -1.0), SVector(0.2, 0.2, 0.5, 1.0, 0.2))
q0 = @isdefined(Q0) ? Q0 : SVector(-0.1, 0.0, -0.5, 0.0)
initial_velocity = @isdefined(INITIAL_VELOCITY) ? INITIAL_VELOCITY : nothing
target_velocity  = @isdefined(TARGET_VELOCITY)  ? TARGET_VELOCITY  : nothing
```

Defines the scenario: the **obstacle** (`_triangle(cx, hw, h)` is the polygon with base
`[cx − hw, cx + hw]` and apex height `h`), the **foothold**, the **discretization**, the **state
box**, the **initial configuration** and the two **velocity ranges** (§3). Each can be overridden
from a driver script.

```julia
X_box = LazySets.Hyperrectangle(; low = x_lb, high = x_ub)
U_box = LazySets.Hyperrectangle(; low = SVector(-disc.u_max, …), high = SVector(disc.u_max, …))
@time removed = RP.infeasible_cells(X_box, disc.state_grid, obstacle, geometry)
margin = RP.swing_foot_deviation_bound(geometry, MP.get_h(disc.state_grid) ./ 2)
X = UT.set_minus(X_box, removed)
```

* `X_box` and `U_box` are the state and input boxes.
* `RP.infeasible_cells` returns the **cells to remove** as an `MP.CellUnion`. The swing foot does not
  depend on the ankle momentum, so it decides on the 4-D configuration cells once and lifts the
  result along `p₁`. A cell is removed if
  * the swing foot of its centre lies below `z = -margin` (it is *surely* below the ground), or
  * the swing foot of its centre is within `margin` of the obstacle (it *might* touch it).
* `margin` is the **Lipschitz margin** `Σᵢ gradᵢ · rᵢ`: the bound on how far the swing foot can move
  within a cell (per-joint gradient bounds `(0, 0.748, 0.576, 0.374, 0.172)` m/rad times the cell
  half-sizes), 7.5 cm by default. It is why the carving is sound, and why the start must be
  farther than this from the obstacle.
* `X = UT.set_minus(X_box, removed)` is the **state domain used by the abstraction**: the box minus
  the removed cells.

### 4.3 Section 3: the abstraction

```julia
@time lgb = LGA.LocalGrowthBound(f_sym, x_sym, u_sym, X, U_box, disc.tstep)
concrete_system = ST.get_system(lgb)
```

`LocalGrowthBound` receives the symbolic dynamics, the state variables, the input variables, the
domain and the sampling time. It differentiates the dynamics symbolically, compiles `f` and `∂f/∂x`
for floats **and** intervals (about 5 s), and builds the Dionysos global growth bound with a
Jacobian bound over the outer box of `X`. `concrete_system` is the **discrete-time system** (RK4
integration over `tstep`) on the carved domain; the control problem is posed on it.

```julia
v0 = initial_velocity === nothing ? (x_lb[1] + x_ub[1]) / 2 : sum(initial_velocity) / 2
x0 = MP.get_coord_by_pos(disc.state_grid, MP.get_pos_by_coord(disc.state_grid, SVector(v0, q0...)))
println("x0 = ", x0, …); println("swing foot at x0: ", RP.swing_foot_position(SVector{4}(x0[2:5]), geometry))
MP.get_pos_by_coord(disc.state_grid, x0) in removed && error("the start cell is carved out …")
```

Builds the **concrete start of the closed-loop run**: `q0` with the momentum at the middle of the
initial range, snapped to the centre of its grid cell. It prints `x0`, the ankle acceleration `ṗ₁`
there and the swing-foot position, then **stops with an error if the start cell was carved out**,
*before* the expensive abstraction.

```julia
extra_attributes = one_joint_per_step ? (; state_input_filter = (x, u) -> count(v -> abs(v) > 1e-9, u) <= 1) : (;)
@time optimizer = LGA.build_abstraction(lgb; state_grid = disc.state_grid, input_grid = disc.input_grid,
    print_level = 1, use_implicit_mapping = true, execution_backend = SY.ThreadedBackend(), extra_attributes...)
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
```

`ONE_JOINT_PER_STEP` becomes an optimizer attribute (a state-input filter that keeps only the inputs
with at most one moving joint). `build_abstraction` then runs the Dionysos
`UniformGridAbstraction` optimizer, with the reach sets of `lgb`:

1. the **state mapping** is created implicitly (cell ids by linearization of the ambient box, low
   memory) on the domain `X`;
2. the **input mapping** is the grid `{-du, 0, du}³` (27 inputs);
3. for **every `(cell, input)` pair**, in parallel threads, `lgb` simulates the cell centre,
   computes the global image, verifies a flow enclosure, bounds the Jacobian locally and returns the
   box image; the cells that box overlaps become the successors, unless the box leaves the domain
   (then the pair has no transition).

It prints the number of states (176 225), inputs (27) and forward images (4 758 075). The returned
`optimizer` **caches the abstraction**. Then

```julia
println("transitions: ", SY.get_n_transitions(abstract_system))
println("local growth bound statistics: ", LGA.stats(lgb))
```

prints the number of transitions and the counters of the local bound: `improved ≈ calls` means the
local bound beat the global one almost everywhere, `enclosure_failures = 0` means the sampling period
is compatible with the dynamics.

### 4.4 Section 4: the footstep problem and the controller

```julia
step_pb = RP.step_problem(concrete_system, disc.state_grid; q0, foothold, geometry, initial_velocity, target_velocity)
```

Builds the **reach problem** (`PR.OptimalControlProblem`): the *initial set* is the box around the
cell of `q0` times the initial velocity range; the *target set* is the union of the cells (found by
inverse kinematics of the swing leg) placing the swing foot at `foothold`, restricted to the target
velocity range. There is no cost: the over-approximation is nondeterministic, so the synthesis
works against the worst case (worst-case number of steps). `step_problem` throws an error if the
start cell is not in the domain.

```julia
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
@time MOI.optimize!(optimizer)
println("footstep success: ", MOI.get(optimizer, MOI.RawOptimizerAttribute("success")))
```

Sets the problem **on the same optimizer** (the abstraction is not recomputed) and solves it. The
solver lifts the initial set to abstract cells (every cell it touches) and the target (only cells
entirely inside), computes the set of cells from which the target is reachable **whatever the
nondeterminism does**, and the controller. `success` is true if the whole initial set is in that
set.

```julia
controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
reached(x) = x ∈ step_pb.target_set
traj = ST.get_closed_loop_trajectory(discrete_time_system, controller, x0, 400; stopping = reached)
```

`controller` is the quantized controller: it maps a concrete state to the cell it lies in and
outputs the input stored for that cell. The **closed loop** simulates the concrete (RK4) system from
`x0` under the controller, for at most 400 steps, stopping when the target is reached or when the
controller has no action at the current state. It is the end-to-end check of the abstraction.

```julia
xs = collect(ST.states(traj)); us = collect(ST.inputs(traj))
println("steps: ", length(xs) - 1, ", reached: ", reached(xs[end]))
…println("ankle along the run: q₁ ∈ ", extrema(x[2] for x in xs), ", p₁ ∈ ", extrema(x[1] for x in xs))
```

Prints the number of steps, whether the target was reached, whether the run stopped *outside the
controller's domain* (no action at the current state), and the ranges of the ankle angle and
momentum along the run.

### 4.5 Section 5: visualization

```julia
if @isdefined(RUN_PLOTS) ? RUN_PLOTS : true
    RP.plot_footstep(traj, disc.state_grid, removed; x_lb, x_ub, tstep = disc.tstep, geometry, obstacle,
                     foothold, file_prefix = joinpath(@__DIR__, "biped_4d_normal_form_$(scenario_name)"),
                     ankle_label = "p₁ (N m s)", title = "…")
end
```

Skipped if `RUN_PLOTS = false`; `RP.plot_footstep` also draws nothing if the controller has no action
at `x0` (`traj` has no input), which is the default outcome. Otherwise it saves, next to the script:

1. a **snapshot** figure (`…_step_normal_form.png`): the robot every fifth step, drawn with the
   ground, the obstacle, the foothold and the forward kinematics of `RP`;
2. the **ankle plot** (`…_step_normal_form_ankle.png`): `q₁` and `p₁` against time;
3. an **animated dashboard** (`…_step_normal_form.gif`, via `DI.animate_trajectory_dashboard`): the
   robot, the state panel `(q₃, q₄)` with the slice of carved cells at the current `(q₁, q₂)` drawn
   as a grey background (the carving does not depend on `p₁`), and the input panel `(v₃, v₄)`.

The same function draws the figures of `biped_4d_underactuated.jl`; only the file prefix, the label of
the ankle coordinate and the title differ.

---

## 5. Diagnosing a failed run

| Symptom | Cause | What to try |
| :--- | :--- | :--- |
| `error: the start cell is carved out …` | the swing foot at `Q0` is within the Lipschitz margin of the obstacle (or below the ground) | move the swing foot farther back (a larger swing-hip angle `q₃`, or a bent knee `q₄ < 0`), refine the grid (the margin scales with the cell size), or set `OBSTACLE = nothing` |
| `error: the cell of the initial configuration q0 = … is not in the state domain` | `Q0` outside `X_BOUNDS` | widen `X_BOUNDS` or change `Q0` |
| `error: no reachable target configuration for foothold …` | no `(q₃, q₄)` inside `X_BOUNDS` puts the foot on the foothold | widen `X_BOUNDS` (`q₃ ≥ …`, `q₄ ≤ …`) or move the foothold |
| `enclosure_failures` large in the statistics | `tstep` too long for the dynamics (`tstep·‖∂f/∂x‖ ≳ 1`): those cells use the global bound | shorter `tstep` |
| `footstep success: false` | no controller can *guarantee* the footstep from the whole initial set | see below |

**When `success` is false**, the usual suspects, from cheapest to most expensive to test:

1. **A narrower initial set**: `INITIAL_VELOCITY = (-0.05, 0.05)` (the controller must then work from
   fewer ankle momenta). The default initial set contains *every* momentum cell of the domain, including
   fast ankle falls that may be unrecoverable.
2. **A finer ankle grid** (`dq1`, `dq̇1`): the ankle is an unstable pendulum and its
   over-approximation inflates every step; with coarse cells the worst-case successors span several
   cells and nothing can be guaranteed.
3. **A different start** (`Q0`) or a **larger domain** (`X_BOUNDS`), for instance more room for the
   ankle to fall.
4. **Check the abstraction**: how large is the set from which the target is reachable?

```julia
ctrl_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
n_ctrl = MP.get_n_state(MP.MappedStateSet(ctrl_set, SY.get_state_mapping(abstract_system)))
n_target = length(step_pb.target_set.positions)
# n_ctrl ≈ n_target means the target is reachable from (almost) nowhere else
```

In an earlier diagnostic with a coarse grid (`dq1 = 0.1`) the controllable set was about the size of
the target itself. This diagnostic has not been repeated for the default settings of this document,
nor has a narrower `INITIAL_VELOCITY` been tried with the obstacle in place. The only synthesis in this
family that reported success is the original example's (`biped_4d_underactuated.jl`, centre
simulation, no obstacle, `INITIAL_VELOCITY = (-0.15, 0.15)`), and there the concrete closed loop
then left the controller's domain after 8 steps without reaching the target: its abstraction is not
certified.
