# Symbolic multibody modelling and local-growth-bound abstraction

This document describes the Julia files of `problems/BipedRobot/4D_model_underactuated/` that are
used to design controllers for the four-link biped:

| File | Module | One-line description |
| :--- | :--- | :--- |
| `Multibody2D.jl` | `Multibody2D` | Symbolic planar multibody kinematics, energies and equations of motion in *normal form*. Generic: nothing in it is specific to the biped. |
| `LocalGrowthAbstraction.jl` | `LocalGrowthAbstraction` | Discrete-time transition systems with **cell-wise (local) growth bounds**, plugged into the Dionysos abstraction pipeline. Generic as well. |
| `robot_problem.jl` | `UnderactuatedRobotProblem` | The biped: **the single definition of the mechanism** (symbolic, with `Multibody2D`), its kinematics and its two 5-D velocity-controlled design models derived from it, and the footstep scenario (carving, target set, grids, plots). **No RigidBodyDynamics.** |
| `BipedFourLink.jl`, `four_link_biped.urdf` | `BipedFourLink` | **Simulation and visualization** of the whole robot with RigidBodyDynamics and MeshCat. Kept separate because it is the only file that needs that dependency; see [`BipedFourLink.md`](BipedFourLink.md). |

The folder contains nothing else besides `Documentation/` (these documents) and `tests/`. Two
modules that used to exist, `BipedModel.jl` and `UnderactuatedLocalAbstraction.jl`, were merged into
`robot_problem.jl` because they duplicated it (§4.6).

The experiment that uses all of them is described in
[`biped_4d_normal_form_experiment.md`](biped_4d_normal_form_experiment.md).

**Contents**

1. [Overview and requirements](#1-overview-and-requirements)
2. [`Multibody2D`](#2-multibody2d)
3. [`LocalGrowthAbstraction`](#3-localgrowthabstraction)
4. [`robot_problem.jl`: the biped, its design models and the footstep scenario](#4-robot_problemjl-the-biped-its-design-models-and-the-footstep-scenario)
5. [Tests and known limitations](#5-tests-and-known-limitations)

---

## 1. Overview and requirements

### 1.1 The pipeline

```
 Multibody2D                              LocalGrowthAbstraction                    Dionysos
 ───────────                              ──────────────────────                    ────────
 frames, points, masses  ──►  normal form  ──►  numbers  ──►  SymbolicDynamics ──►  LocalGrowthBound
 (symbolic, in q)             x = (q,v_a,p_u)   substituted    (compiled f, ∂f/∂x      (over-approximation
                              ẋ = f(x) + g(x)τ  late           for floats & intervals)  of every cell image)
                                                                                          │
                                                                                          ▼
                                                              build_abstraction ──► UniformGridAbstraction
                                                                                    optimizer (cached transition
                                                                                    system) ──► control problem
```

The point of the design is that the model stays **symbolic** until the last moment (no
RigidBodyDynamics, whose symbolic handling is slow and unmaintained), so a general robot can be
described once and its parameters (lengths, masses, inertias, gravity) replaced by numbers
afterwards. The abstraction is a **certified over-approximation** (up to the boundary
convention of §3.4) that improves on the global growth bound of Dionysos wherever the dynamics
vary over the state space.

### 1.2 Requirements and how to load the modules

The modules are plain files, not a package. Load them with `include` and `import`/`using`:

```julia
using Symbolics, IntervalArithmetic      # symbolic model; the interval extension of Dionysos
import Dionysos, LazySets, MathOptInterface as MOI
using StaticArrays

dir = joinpath(dirname(dirname(pathof(Dionysos))), "problems", "BipedRobot", "4D_model_underactuated")
include(joinpath(dir, "robot_problem.jl"));          import .UnderactuatedRobotProblem as RP   # the biped; includes Multibody2D.jl
include(joinpath(dir, "LocalGrowthAbstraction.jl")); import .LocalGrowthAbstraction as LGA
const MB = RP.Multibody2D                # the symbolic multibody module, to build anything else
```

`robot_problem.jl` includes its own copy of `Multibody2D.jl` (reachable as `RP.Multibody2D`), so
objects of one copy must not be passed to functions of another: when you work with the biped, use
`RP.Multibody2D` and do not include `Multibody2D.jl` a second time. To model another robot, include
`Multibody2D.jl` on its own (`include(joinpath(dir, "Multibody2D.jl")); using .Multibody2D`).
`LocalGrowthAbstraction.jl` does not depend on `Multibody2D.jl`: it only needs symbolic
expressions of the dynamics.

Packages needed in the active environment: `Dionysos`, `Symbolics`, `IntervalArithmetic`,
`MathematicalSystems`, `LazySets`, `StaticArrays`, `MathOptInterface`, and `Plots` for
`robot_problem.jl`. **RigidBodyDynamics is not needed** by anything except `BipedFourLink.jl`. The
code was developed and tested with Julia 1.12, Symbolics 7.45 and IntervalArithmetic 1.0. The
`problems/BipedRobot` project lists all of them (with RigidBodyDynamics and MeshCat, needed by
`BipedFourLink.jl` only).

### 1.3 Conventions used throughout

* **Angles are counter-clockwise** (radians), with `x` pointing forwards and `z`/`y` up. A
  rotation by `θ` is `R(θ) = [cos θ -sin θ; sin θ cos θ]`.
* Every symbolic object is a `Symbolics.Num`. Vectors are `Vector{Num}`, matrices `Matrix{Num}`.
* Dionysos conventions are followed: sets are `LazySets` objects (a `Hyperrectangle` is built
  from bounds as `Hyperrectangle(; low, high)`), grids are `MP.GridFree(origin, step)`, and
  every function that returns a sequence of coordinates uses `SVector`.

---

## 2. `Multibody2D`

Symbolic planar multibody kinematics, energies and equations of motion.

### 2.1 Concepts

* A **`MultibodySystem`** holds the generalized coordinates `q`, their velocities `q̇` (vectors of
  symbolic variables) and the gravity vector.
* A **`Frame`** is a 2D reference frame given by a translation and an angle with respect to a
  parent frame (the world frame or another frame). Frames form a tree.
* A **`Point`** is a location in a frame, optionally carrying a rigid body (mass, centre of mass,
  moment of inertia about the centre of mass).
* Translations and angles may be numbers, symbolic parameters, or expressions of `q`.
* **Everything time-varying must be written in terms of `q`.** There is no time variable: the
  velocity of any quantity is `∂(·)/∂q · q̇`.

### 2.2 Running example

All the examples of this section use the double pendulum below (link `i` has length `lᵢ`, its
centre of mass at `lcᵢ`, mass `mᵢ`, inertia `Iᵢ`).

```julia
using Symbolics
include("Multibody2D.jl"); using .Multibody2D

@variables q[1:2] qd[1:2] l1 l2 lc1 lc2 m1 m2 I1 I2 g
sys   = MultibodySystem(collect(q), collect(qd); gravity = g)
link1 = add_frame!(sys, :link1, sys.world, [0, 0], q[1])
link2 = add_frame!(sys, :link2, link1,     [l1, 0], q[2])
add_point!(sys, :body1, link1, [lc1, 0]; mass = m1, inertia = I1)
add_point!(sys, :body2, link2, [lc2, 0]; mass = m2, inertia = I2)
tip = add_point!(sys, :tip, link2, [l2, 0])
```

### 2.3 Reference

#### `MultibodySystem(q, q̇; gravity = 9.81)`

*Description.* The container of the model: `q` and `q̇` are equal-length vectors of symbolic
variables (`collect(q)` for `@variables q[1:n]`). `gravity` is either a scalar `g` (acceleration
`[0, -g]`) or a 2-vector in world coordinates; it may be symbolic. Fields: `q`, `qd`, `gravity`,
`world` (the root `Frame`), `frames`, `points` (both in insertion order).

*Use case.* Start every model with it. `sys[:name]` returns the frame or point with that name
(`sys[:world]` is the world frame); names must be unique.

```julia
sys                       # MultibodySystem(2 coordinates, 2 frames, 3 points)
sys[:tip] === tip         # true
MultibodySystem([q[1]], [qd[1], qd[2]])   # throws DimensionMismatch (different lengths)
```

#### `Frame` and `add_frame!(sys, name, parent, translation, angle)`

*Description.* `Frame` is the type of a reference frame. A point with coordinates `x` in the frame
has coordinates `translation + R(angle) * x` in its parent. Always create frames with
`add_frame!`, which checks that the name is unused and that `parent` belongs to `sys`.
`translation` is a 2-vector in the coordinates of `parent`; `angle` is the rotation relative to
`parent`. Both may be symbolic and may depend on `q`.

*Use case.* A kinematic chain (each frame hangs from the previous one), or a frame whose origin
moves with a generalized coordinate, e.g. a cart:

```julia
cart = add_frame!(sys, :cart, sys.world, [x, 0], 0)      # translation depends on the coordinate x
pole = add_frame!(sys, :pole, cart, [0, 0], θ)
```

#### `Point` and `add_point!(sys, name, frame, position; mass, com, inertia)`

*Description.* A point at `position` (2-vector, coordinates in `frame`). Without `mass` it is a
purely kinematic marker (a foot, a joint location). With `mass` a rigid body is attached to
`frame`: `com` is the position of its centre of mass **in the coordinates of `frame`** (default:
the point itself, i.e. a point mass) and `inertia` its moment of inertia about the centre of mass
(default `0`). `com` and `inertia` without `mass` throw an `ArgumentError`.

*Use case.* Links are rods whose centre of mass is not at the frame origin; motors are point
masses with zero inertia; a foot is a marker:

```julia
add_point!(sys, :rod,   link1, [0, 0];  mass = m1, com = [lc1, 0], inertia = m1*l1^2/12)
add_point!(sys, :motor, link2, [0, 0];  mass = mM)               # point mass at the frame origin
foot = add_point!(sys, :foot, link2, [l2, 0])                    # marker, no mass
```

#### `relative_transform(from, to) -> (t, a)`

*Description.* The planar transform from the coordinates of frame `from` to those of frame `to`:
a point `x` in `from` is `t + R(a) * x` in `to`. The path goes through the lowest common ancestor
of the two frames, so `a` is expressed with the fewest terms. The returned expressions are *not*
simplified (you may see `cos(0)*l1`); apply `trig_reduce` if you need them clean.

*Use case.* Where is link 2 seen from link 1, e.g. to place a sensor.

```julia
t, a = relative_transform(link2, link1)     # t = [l1, 0] (unsimplified), a = q[2]
```

#### `rotation_matrix(θ)`

*Description.* The 2×2 `Matrix{Num}` `R(θ)`.

*Use case.* Rotate a vector between frames by hand, e.g. the world velocity into a body frame.

```julia
rotation_matrix(q[1])     # [cos(q[1]) -sin(q[1]); sin(q[1]) cos(q[1])]
```

#### `position(p, frame = root frame; simplify = true)`

*Description.* Coordinates of the point `p` in `frame` (default: the world frame). It extends
`Base.position`. `simplify = true` applies `trig_reduce`.

*Use case.* Forward kinematics, e.g. to constrain the height of a foot or to plot.

```julia
position(tip)           # [cos(q[1])*l1 + cos(q[1]+q[2])*l2,  l1*sin(q[1]) + l2*sin(q[1]+q[2])]
position(tip, link1)    # [l1 + cos(q[2])*l2,  l2*sin(q[2])]
```

#### `velocity(sys, p; relative_to = sys.world, expressed_in = relative_to, simplify = true)`

*Description.* The velocity `∂(position)/∂q · q̇` of `p`, **as seen by an observer attached to
`relative_to`**, with components in the axes of `expressed_in`. The two keywords are different:
`relative_to` selects who observes (a point at rest in that frame has zero velocity),
`expressed_in` only selects the axes in which the vector is written.

*Use case.* World velocity of a foot (e.g. to impose a no-slip condition), or the velocity of a
joint relative to its parent link.

```julia
velocity(sys, tip)                          # world velocity, world axes
velocity(sys, tip; relative_to = link1)     # [-l2*qd[2]*sin(q[2]),  l2*qd[2]*cos(q[2])]
velocity(sys, tip; expressed_in = link2)    # world velocity written in the axes of link2
```

#### `mass_matrix(sys; simplify = true) -> Matrix{Num}`

*Description.* The symmetric matrix `D(q)` such that the kinetic energy is `T = ½ q̇ᵀ D(q) q̇`.
Each body contributes `m JᵀJ + I J_θᵀ J_θ`, with `J = ∂(centre of mass)/∂q` and
`J_θ = ∂(frame angle)/∂q`, both in the world frame.

*Use case.* Any Lagrangian-mechanics computation, or a check against another model.

```julia
D = mass_matrix(sys)
D[1, 2]        # I2 + lc2^2*m2 + cos(q[2])*l1*lc2*m2     (the textbook double-pendulum entry)
```

#### `kinetic_energy(sys; simplify = true)`, `potential_energy(sys; simplify = true)`, `lagrangian(sys; simplify = true)`

*Description.* `T = ½ q̇ᵀ D q̇`; `V = -Σ m gᵀ r_com` over all bodies (gravity as given to
`MultibodySystem`); `L = T - V`. All are scalar `Num`s.

*Use case.* Energy-based checks (conservation along a simulation), or the starting point of other
derivations (`Symbolics.derivative(lagrangian(sys), q[1])` is `∂L/∂q₁`).

```julia
potential_energy(sys)     # (g*l1*m2 + g*lc1*m1)*sin(q[1]) + g*lc2*m2*sin(q[1]+q[2])
```

#### `bias_vector(sys; simplify = true) -> Vector{Num}`

*Description.* `h(q, q̇) = ∂L/∂q - Ḋ(q) q̇` (Coriolis, centrifugal and gravity terms, sign included),
such that the Euler–Lagrange equations read `D(q) q̈ = h(q, q̇) + B τ`.

*Use case.* Second-order dynamics `q̈ = D⁻¹ (h + Bτ)` for a simulator, or to validate the
normal form (which is built from the same ingredients).

```julia
bias_vector(sys)[2]       # -cos(q[1]+q[2])*g*lc2*m2 - l1*lc2*m2*qd[1]^2*sin(q[2])
```

#### `trig_reduce(x)`

*Description.* Rewrites a symbolic expression (or each entry of an array) as a *linear*
combination of `sin(arg)`/`cos(arg)` with canonical arguments, using the product-to-sum
identities. This is what turns `cos(q1)cos(q1+q2) + sin(q1)sin(q1+q2)` into `cos(q2)` and makes
every `sin² + cos²` cancel, which `Symbolics.simplify` alone does not do. Non-trigonometric
structure (parameters, polynomials in `q̇`) is only expanded. All the functions above call it by
default (`simplify = true`).

*Use case.* Compare two expressions that are equal only through trigonometric identities, or
clean up an expression you built by hand:

```julia
trig_reduce(cos(q[1])*cos(q[1]+q[2]) + sin(q[1])*sin(q[1]+q[2]))     # cos(q[2])
Symbolics.simplify(cos(q[1])*cos(q[1]+q[2]) + sin(q[1])*sin(q[1]+q[2]))   # unchanged
iszero(trig_reduce(lagrangian(sys) - kinetic_energy(sys) + potential_energy(sys)))   # true
```

Pass `simplify = false` to any of the functions to skip it (faster, uglier expressions).

#### `selection_matrix(n, actuated) -> Matrix{Num}`

*Description.* The `n × length(actuated)` input matrix `B` in which input `j` acts directly on
coordinate `actuated[j]`: column `j` is the `actuated[j]`-th Euclidean basis vector, so the rows of
the unactuated coordinates are zero. Repeated or out-of-range indices throw an `ArgumentError`.

*Use case.* The biped has four coordinates and torque only on the last three:

```julia
B = selection_matrix(4, [2, 3, 4])        # [0 0 0; 1 0 0; 0 1 0; 0 0 1]
```

#### `NormalForm` and `normal_form(sys, B; input = nothing, momenta = nothing, simplify = true)`

*Description.* The dynamics of `sys` in **normal form**. `B` is the `n × m` input matrix of the
Euler–Lagrange equations

```
d/dt(∂L/∂q̇ᵢ) − ∂L/∂qᵢ = Bᵢ τ ,      L = T − V .
```

A coordinate whose row of `B` is zero is *unactuated*; the others are *actuated* (index sets `a`
and `u`). The state is `x = (q, v_a, p_u)` with `v_a = q̇_a` the actuated velocities and
`p_u = (∂L/∂q̇)_u = (D q̇)_u` the **generalized momenta of the unactuated coordinates**, and

```
q̇_a = v_a ,      q̇_u = D_uu⁻¹ (p_u − D_ua v_a) ,
ṗ_u = ∂L/∂q_u ,                                          (independent of τ, since B_u = 0)
S v̇_a = h_a − D_au D_uu⁻¹ h_u + B_a τ ,    S = D_aa − D_au D_uu⁻¹ D_ua   (Schur complement),
```

with `h = ∂L/∂q − Ḋ q̇` evaluated at the `q̇` above. For one or two unactuated coordinates `S` is
written as `N/den` with a trigonometric-polynomial numerator, so `S⁻¹` has no nested fractions.

The result is a `NormalForm` with fields

| Field | Content |
| :--- | :--- |
| `state` | `x = [q; v_a; p_u]` (the entries of `v_a` are `sys.qd[actuated]`) |
| `input` | the inputs `τ` (default variables `τ₁, τ₂, …`) |
| `actuated`, `unactuated` | the two index sets |
| `velocities` | the full `q̇` as a function of the state |
| `xdot` | `ẋ = drift + gain * input` |
| `drift`, `gain` | the `τ`-independent part of `ẋ` and its input matrix |

`input` and `momenta` optionally supply your own symbolic variables for `τ` and `p_u` (default
momenta: `p_i`, with `i` the coordinate index).

*Use case.* The acrobot (double pendulum actuated at the elbow only). The momentum of the
unactuated shoulder replaces its velocity, and its dynamics contain no torque:

```julia
nf = normal_form(sys, selection_matrix(2, [2]))
nf                 # NormalForm(4 states, 1 inputs; actuated = [2], unactuated = [1])
nf.state           # [q[1], q[2], qd[2], p₁]
nf.xdot[end]       # ṗ₁ = cos(q[1])*(-g*l1*m2 - g*lc1*m1) - cos(q[1]+q[2])*g*lc2*m2     (no τ₁)
nf.gain[3]         # (D₁₁) / (I1*I2 + … − ½cos(2q[2])*l1²*lc2²*m2²) : the gain of τ on q̈₂, = D₁₁/det D
```

Systems with one, two or three unactuated coordinates and fully actuated systems (`B` with no zero
row) are covered by the tests; a system with *no* actuated coordinate is not.

##### From the symbolic model to numbers

Because everything is a `Num`, parameters can be substituted after the model is built. Compile the
result with **`cse = true`** (common-subexpression elimination); without it a three-actuator
normal form takes minutes to compile (a 4-link chain: 178 s to the first call, against 0.07 s):

```julia
vals = Dict(l1 => 1.0, l2 => 1.0, lc1 => 0.5, lc2 => 0.5, m1 => 1.0, m2 => 1.0,
            I1 => 0.083, I2 => 0.083, g => 9.81)
xdot = Symbolics.substitute.(nf.xdot, Ref(vals))             # still symbolic in x and τ
f, f! = build_function(xdot, nf.state, nf.input; expression = Val(false), cse = true)
f([0.1, 0.2, 0.3, 0.0], [0.5])                               # ẋ at x = (q₁, q₂, q̇₂, p₁), τ = 0.5
```

Alternatively keep the parameters as run-time arguments:
`build_function(nf.xdot, nf.state, nf.input, params; expression = Val(false), cse = true)`.

### 2.4 Notes and limitations

* Planar systems only; the generalized coordinates are your choice (one per degree of freedom).
* `normal_form` solves the linear systems symbolically by Gaussian elimination *without pivoting*,
  which is valid because every block it inverts is symmetric positive definite. With three
  actuated coordinates and fully symbolic parameters the printed expressions are huge
  (≈ 1.4 M characters for a 4-link chain); that is why numbers are substituted and `cse = true`
  is used.
* The kinematic functions are exact and cheap; `normal_form` of a 4-DoF chain with symbolic
  parameters takes about 1 s once compiled.

---

## 3. `LocalGrowthAbstraction`

Discrete-time transition systems with **cell-wise (local) growth bounds**.

### 3.1 Why

The classical growth bound of Dionysos (`ST.ContinuousTimeGrowthBound`) uses one Jacobian bound
`L(u)` valid on the whole state space `X`: a cell of radius `r` is mapped to `Φ(c, u) ± ρ`, with
`ρ̇ = L(u) ρ`, `ρ(0) = r`. If the dynamics are strongly non-linear only in a small part of `X`,
`L` is dictated by the worst point and the image of every other cell is needlessly large, which
creates **spurious transitions**.

### 3.2 The algorithm (for one cell `R = c ± r` and input `u`, sampling time `τ`)

1. **Centre.** Simulate the centre with the RK4 integrator of Dionysos: `Φ(c, u)`.
2. **Global image.** `Φ(c, u) ± ρ_glob`, where `ρ_glob` comes from the Dionysos global growth
   bound, whose Jacobian bound is taken with interval arithmetic over all of `X`.
3. **Larger cell, local bound.** Take a box `B ⊇ R ∪ image`, **verified** to enclose the whole
   flow tube (see below), and bound the Jacobian over `B` only, again with interval arithmetic:
   the local growth bound `L_B`.
4. **Local image.** Integrate `ρ̇ = L_B ρ` and take `Φ(c, u) ± min(ρ_glob, ρ_loc)`.

*Why `B` is verified.* The growth-bound theorem needs `L_B` to bound the Jacobian on a convex set
containing the trajectories of **all** points of `R` during the **whole** sampling period, not
just at its end points. `B` is therefore checked with an interval Picard–Lindelöf test,
`R + [0, τ]·f(B, u) ⊆ B`, and enlarged until it passes. If no such box is found (this happens
when `τ·‖∂f/∂x‖` is not small, roughly `≳ 1`) the cell falls back to the global bound, which is
always valid. Both images are valid, so their intersection is valid, and the local image is never
worse than the global one. Optionally steps 3–4 are repeated (`refinements`), each time starting
from the previous, smaller image.

The Jacobian bound is the *Metzler* majorant required by the radius ODE: the supremum of the
diagonal entries and the magnitude of the off-diagonal ones.

### 3.3 Reference

#### `SymbolicDynamics(fsym, x, u)` and `SymbolicDynamics(system)`

*Description.* The dynamics `ẋ = f(x, u)` given symbolically and compiled once into three
callables (fields):

| Field | Signature | Use |
| :--- | :--- | :--- |
| `f` | `(x, u) -> SVector{N}` | plain floating-point evaluation (centre simulation) |
| `f_interval` | `(xbox, ubox) -> Vector` | interval enclosure of `f` over boxes (flow enclosure) |
| `jacobian_interval` | `(xbox, ubox) -> Matrix` | interval enclosure of `∂f/∂x` (growth bounds) |

`fsym` is a vector of `Symbolics` expressions in the variable vectors `x` and `u` (for instance the
numerical normal form of `Multibody2D`). The second constructor takes a Dionysos
`ConstrainedBlackBoxControlContinuousSystem` and *traces* its `f` with symbolic variables, so `f`
must be written with generic operations. Compilation uses `nanmath = false` (the default
NaNMath-backed `sin`/`cos` do not accept intervals) and `cse = true`.

*Use case.* Bring any smooth dynamics into the interval machinery, e.g. a pendulum:

```julia
@variables x[1:2] u[1:1]
fp  = [x[2], -9.81*sin(x[1]) - 0.5*x[2] + u[1]]
dyn = SymbolicDynamics(fp, collect(x), collect(u))
dyn.f([0.1, 0.2], [0.0])            # [0.2, -1.0794]
```

#### `jacobian_bound(dynamics, xlow, xhigh, u) -> SMatrix`

*Description.* The Metzler growth-bound matrix `L` of `∂f/∂x`, bounded with interval arithmetic
for all states in the box `[xlow, xhigh]` and the input `u`. Entries may be `Inf` if the interval
evaluation is unbounded (e.g. a denominator whose interval contains zero).

*Use case.* Inspect how much a local box tightens the bound. For the pendulum, `|∂(−g sin x₁)/∂x₁|`
is `g` over a full turn but small near `x₁ = π/2`:

```julia
jacobian_bound(dyn, [-3.0, -6.0], [3.0, 6.0], [0.0])   # [0 1; 9.81 -0.5]   (whole state space)
jacobian_bound(dyn, [1.4, -1.0], [1.7, 1.0], [0.0])    # [0 1; 1.67 -0.5]   (a small box)
```

#### `flow_enclosure(dynamics, cell_low, cell_high, guess_low, guess_high, u, tstep; max_iterations = 10, inflation = 0.1, shrink = false)`

*Description.* Returns `(lo, hi)`, a box that provably contains the trajectory `x(t)`,
`t ∈ [0, tstep]`, of every initial state in the cell `[cell_low, cell_high]` under the constant
input `u`, or `nothing` if none was found. The search starts from the hull of the cell and
`guess` (the conservative image); when the Picard test fails the box is replaced by the hull of
itself and `R + [0, τ]·f(B, u)`, inflated by `inflation` times its width, up to `max_iterations`
times. With `shrink = true` the tighter (still valid) `R + [0, τ]·f(B, u)` is returned. The
endpoint arithmetic is plain floating point, not outward rounded.

*Use case.* Check that a sampling period is compatible with the dynamics before building an
abstraction (`nothing` means the local bound will fall back to the global one):

```julia
cell = LazySets.Hyperrectangle(SVector(1.5, 0.0), SVector(0.1, 0.2))
glob = ST.get_over_approximation_map(lgb.global_growth)(cell, SVector(0.0))
flow_enclosure(dyn, LazySets.low(cell), LazySets.high(cell),
               LazySets.low(glob), LazySets.high(glob), [0.0], 0.1)    # (lo, hi) or nothing
```

#### `LocalGrowthBound` (a `Dionysos.System.DiscreteTimeSystemOverApproximation`)

*Description.* The over-approximation of the sampled dynamics (sampling period `tstep`) with local
growth bounds. Three constructors:

```julia
LocalGrowthBound(fsym, x, u, X, U, tstep; kwargs...)        # from symbolic dynamics and the sets
LocalGrowthBound(system, tstep; kwargs...)                  # from a Dionysos continuous system (traced)
LocalGrowthBound(system, dynamics::SymbolicDynamics, tstep; kwargs...)
```

`X` and `U` are `LazySets` sets (`X` may be a box minus carved cells: the global Jacobian bound uses
its outer box). The cell centre is simulated with the `f` of the system (RK4, via `ST.discretize`).
Keywords:

| Keyword | Default | Meaning |
| :--- | :--- | :--- |
| `refinements` | `1` | how many times steps 3–4 are applied (`0` = the global growth bound of Dionysos) |
| `ngrowthbound` | `ST.DEFAULT_NUM_SUBSTEPS` (5) | RK4 substeps for the radius ODE `ρ̇ = Lρ` |
| `nsystem` | `ST.DEFAULT_NUM_SUBSTEPS` (5) | RK4 substeps for the centre simulation |
| `max_enclosure_iterations` | `10` | iterations of the Picard enclosure before falling back |
| `enclosure_inflation` | `0.1` | relative inflation of the box at every failed iteration |
| `shrink_enclosure` | `false` | use the tighter verified box `R + [0,τ]·f(B,u)` for the Jacobian |
| `boundary_tolerance` | `1e-9` | see §3.4 |

Fields of interest: `system` (the discrete-time RK4 system on `X`, what a control problem
refers to), `global_growth` (the Dionysos `DiscreteTimeGrowthBound`, also the fallback and the
baseline for comparisons), `dynamics`, `tstep`, `stats`.

*Use case.* The pendulum on a grid, compared with the plain global bound:

```julia
X = LazySets.Hyperrectangle(; low = [-3.0, -6.0], high = [3.0, 6.0])
U = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
lgb = LocalGrowthBound(fp, collect(x), collect(u), X, U, 0.1)

cell = LazySets.Hyperrectangle(SVector(1.5, 0.0), SVector(0.1, 0.2))
over_local  = ST.get_over_approximation_map(lgb)(cell, SVector(0.0))
over_global = ST.get_over_approximation_map(lgb.global_growth)(cell, SVector(0.0))
LazySets.radius_hyperrectangle(over_global)    # [0.1247, 0.2971]
LazySets.radius_hyperrectangle(over_local)     # [0.1212, 0.2245]
```

As a Dionysos object it implements `ST.get_system`, `ST.input_cache` (the global radius is
integrated once per input and reused for every cell), `ST.reach_set` and
`ST.get_over_approximation_map`, so it can be used wherever Dionysos expects an
over-approximation.

#### `LocalGrowthBoundStats`, `stats(lgb)`, `reset_stats!(lgb)`

*Description.* Thread-safe counters of a `LocalGrowthBound`: `calls` (cell/input pairs processed),
`enclosure_failures` (no verified enclosure, global bound used), `unbounded` (the local Jacobian
bound was not finite, global bound used) and `improved` (pairs whose local radius is strictly
smaller than the global one in at least one dimension). `stats(lgb)` returns a `NamedTuple`;
`reset_stats!(lgb)` zeroes them.

*Use case.* Diagnose an abstraction: many `enclosure_failures` mean the sampling period is too
long for the dynamics, `improved ≈ calls` means the local bound helps everywhere.

```julia
stats(lgb)     # (calls = 2523, enclosure_failures = 0, unbounded = 0, improved = 2172)
reset_stats!(lgb)
```

#### `build_abstraction(approx; state_grid, input_grid, print_level = 0, abstraction_file = nothing, abstraction_info = nothing, attributes...) -> optimizer`

*Description.* Builds the discrete-time transition system (symbolic model) of the system of
`approx` with the Dionysos `UniformGridAbstraction` optimizer in `USER_DEFINED` mode, using the
reach sets of `approx`. Any `DiscreteTimeSystemOverApproximation` works: a `LocalGrowthBound`, or
`lgb.global_growth` as the baseline. `state_grid` and `input_grid` are Dionysos grids; every extra
keyword is forwarded as an optimizer attribute (`execution_backend`, `use_implicit_mapping`,
`state_input_filter`, …). It returns the optimizer, which **caches the abstraction**: a control
problem can be set on it afterwards. With `abstraction_file`, the transitions are cached on disk: loaded
if the file exists (`load_abstraction!`, which refuses a file whose `abstraction_info`, number of states
and inputs, or state numbering differ), otherwise computed and written (`save_abstraction`: JLD2, three
integer vectors, 12 bytes per transition). This is how an abstraction computed on a server (see
`slurm_array/README.md`) is used elsewhere; the compiled dynamics are not saved, only the transitions the
control synthesis needs. With a `SY.SlurmArrayBackend` the call computes one chunk of the transitions
of a Slurm array, and `merge_chunks = (; outdir, nchunks)` reads the chunks instead of computing
(`merge_chunks!`).

*Use case.* Build the abstraction, inspect it, then synthesize a controller on the same optimizer:

```julia
const MP = Dionysos.Mapping; const SY = Dionysos.Symbolic
opt = build_abstraction(lgb;
        state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.2, 0.4)),
        input_grid = MP.GridFree(SVector(0.0), SVector(2.0)),
        execution_backend = SY.ThreadedBackend())          # forwarded to the optimizer
abs_sys = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
SY.get_n_state(abs_sys), SY.get_n_input(abs_sys), SY.get_n_transitions(abs_sys)   # (841, 3, 11573)

# the same grids with the global bound: more (spurious) transitions
base = build_abstraction(lgb.global_growth; state_grid = …, input_grid = …)      # 12211 transitions

# control synthesis on the cached abstraction (a reach problem; set the problem, then optimize)
PR = Dionysos.Problem
problem = PR.OptimalControlProblem(ST.get_system(lgb), initial_box, target_box,
                                   nothing, nothing, PR.Infinity())
MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
MOI.optimize!(opt)
controller = MOI.get(opt, MOI.RawOptimizerAttribute("concrete_controller"))
```

### 3.4 The `boundary_tolerance` convention

If a coordinate evolves **exactly on the grid** (for instance `q̇ = u` with `τ·u` a multiple of the
cell size), the image of a cell is exactly a cell and its boundary touches the neighbours along a
face. Dionysos treats cells as closed (its `OUTER` lookup has `tol = 0`), so every such face
contact would count as an overlap and each motor joint would get three successors instead of one,
which makes controller synthesis hopeless. `LocalGrowthBound` therefore deflates the radius of
every image by `boundary_tolerance × (cell radius)` (default `1e-9`): the half-open-cell
convention of grid abstractions. It is far above floating-point noise (≈ 1e-17) and far below any
real growth, but it means that an image protruding into a neighbour by less than that is ignored.
Use `boundary_tolerance = 0` for the strict closed-cell behaviour.

### 3.5 Notes and limitations

* Continuous-time, noise-free systems only.
* The enclosure test needs `τ·‖∂f/∂x‖` not large. In the biped (`‖∂f/∂x‖ ≈ 9`) a sampling period
  of 0.05 s works and 0.1 s is marginal; with high-gain dynamics every cell falls back to the
  global bound (visible in `stats`).
* The endpoint arithmetic of the enclosure is not outward rounded, and `RK4` simulates the centre;
  as everywhere in Dionysos, "certified" means up to numerical integration and rounding error.
* Cost: about 0.1 ms of interval arithmetic per cell/input pair (the biped: 4.8 M pairs in
  ≈ 2 minutes on 10 threads). It is thread safe: use `SY.ThreadedBackend()`.

---

## 4. `robot_problem.jl`: the biped, its design models and the footstep scenario

Module `UnderactuatedRobotProblem` (`RP` below). The biped is defined **once**, symbolically, with
`Multibody2D`; everything else in the file is derived from that definition (the kinematics, the
normal-form design model) or is the footstep scenario built on top of it (obstacle carving, target set,
grids, start state, closed loop, plots). There is no second copy of the mechanism to keep consistent by
hand. It needs Symbolics, Dionysos, LazySets, MathOptInterface, StaticArrays and Plots: **no
RigidBodyDynamics**. `examples/BipedRobot/biped_4d_normal_form.jl` builds on it.

```
 symbolic_model()  ──►  Geometry(model, values)    forward kinematics (compiled) ──► carving, target set, plots
 (Multibody2D,     ──►  reduced_dynamics ─► substitute_parameters ─► x = (p₁, q), symbolic f ──► local-growth-bound abstraction
  symbolic
  parameters)
```

The biped (ankle `q₁` unactuated, grounded knee `q₂`, hip `q₃` and swing knee `q₄` actuated):

```
 world origin = ankle (RF1, q₁, unactuated)
   └─ tibia, length L2 ─► grounded knee (RF2, q₂)
        └─ femur, length L1 ─► hip (RF3, q₃)
             └─ swing femur, length L1 (hanging down) ─► swing knee (RF4, q₄)
                  └─ swing tibia, length L2 ─► swing foot
```

### 4.1 The mechanism, and its consistency with the URDF

The links are **uniform thin rods** (moment of inertia `m L²/12` about their centre of mass, which is
at the middle of the link) and the motors are **point masses** with no rotor inertia; the masses sit
at the joints. `four_link_biped.urdf` describes exactly this robot, so the symbolic design model and
the RigidBodyDynamics simulation of `BipedFourLink.jl` are the same mechanism, which is checked three
ways: the `q̈₁` of the velocity-controlled model (written with the ankle velocity, `q̈₁ = h₁/D₁₁`, in the
tests) agrees with RigidBodyDynamics on the URDF to **machine precision** (reference values kept in
`tests/test_robot_problem.jl`, so the guarantee survives the removal of the dependency); the torque-level accelerations `D q̈ = h + Bτ` agree; and the
forward kinematics of the URDF (joint and foot positions) agree with `Geometry`
(`tests/test_biped_four_link.jl`).

### 4.2 The mechanism and its numbers

#### `symbolic_model() -> NamedTuple`

*Description.* Builds the biped as a `MultibodySystem` with **symbolic** parameters and its normal
form. Returns `(; sys, nf, B, parameters)` (`parameters` is the `NamedTuple` of the symbolic
variables `(L1, L2, m1, m2, mM, g)`). The frames (all translations in the parent frame, all angles
relative to the parent, counter-clockwise) are

| Frame | Parent | Translation | Angle | Meaning |
| :--- | :--- | :--- | :--- | :--- |
| `RF1` | world | `[0, 0]` | `q₁` | grounded tibia, ankle at the origin |
| `RF2` | `RF1` | `[0, L2]` | `q₂` | grounded knee |
| `RF3` | `RF2` | `[0, L1]` | `q₃` | hip |
| `RF4` | `RF3` | `[0, -L1]` | `q₄` | swing knee |

and the points are

| Point | Frame | Position (= centre of mass) | Mass | Inertia |
| :--- | :--- | :--- | :--- | :--- |
| `link1` | `RF1` | `[0, L2/2]` | `m2` | `m2 L2²/12` |
| `link2` | `RF2` | `[0, L1/2]` | `m1` | `m1 L1²/12` |
| `link3` | `RF3` | `[0, -L1/2]` | `m1` | `m1 L1²/12` |
| `link4` | `RF4` | `[0, -L2/2]` | `m2` | `m2 L2²/12` |
| `grounded_knee` | `RF2` | `[0, 0]` | `mM` | `0` |
| `hip` | `RF3` | `[0, 0]` | `2 mM` | `0` |
| `swing_knee` | `RF4` | `[0, 0]` | `mM` | `0` |
| `swing_foot` | `RF4` | `[0, -L2]` | `0` | `0` |

The input matrix is `B = [0 0 0; I₃]` (the last three coordinates are actuated). Building it takes a
few seconds, most of it the normal form: build it **once** and pass it to the functions below.

*Use case.* A model you can inspect or reuse symbolically (kinematics of the swing foot, energies):

```julia
model = RP.symbolic_model()
model.nf                                                   # NormalForm(8 states, 3 inputs; actuated = [2, 3, 4], unactuated = [1])
RP.Multibody2D.position(model.sys[:swing_foot])            # swing-foot position in the world frame, symbolic
```

#### `default_parameters() -> NamedTuple`

*Description.* Numerical values of the parameters, those of `four_link_biped.urdf`: `L1 = 0.202` m
(femur), `L2 = 0.172` m (tibia), `m1 = 0.0202` kg (femur), `m2 = 0.0172` kg (tibia), `mM = 0.4` kg
(motors, by far the dominant mass; the hip carries `2 mM`), `g = 9.81`.

*Use case.* A baseline you can modify field by field:

```julia
vals = merge(RP.default_parameters(), (; mM = 0.6, L2 = 0.18))   # a heavier motor, a longer tibia
```

#### `substitute_parameters(model, expressions, values = default_parameters())`

*Description.* Replaces the symbolic parameters of `model` by the numbers `values` (a `NamedTuple`
with the same keys) in a symbolic expression or an array of them. What is left is symbolic in the
coordinates (and inputs) only. It is the single place where numbers enter the chain: the step
"build the model symbolically, put the numbers in later" of §2.3, used by `Geometry`,
the kinematics and by the abstraction.

*Use case.*

```julia
D11 = RP.substitute_parameters(model, RP.Multibody2D.mass_matrix(model.sys)[1, 1], vals)   # D₁₁(q), numbers in
```

### 4.3 Kinematics

#### `Geometry(model = symbolic_model(), values = default_parameters())`, `default_geometry()`

*Description.* The two segment lengths (`tibia = values.L2`, `femur = values.L1`) and the forward
kinematics of the joints, **compiled from the symbolic positions** of `model` (grounded knee, hip,
swing knee, swing foot): the carving, the target set and the plots cannot disagree with the dynamics.
`default_geometry()` is `Geometry()`, i.e. the numbers of `default_parameters()` (it builds its own
symbolic model: pass the one you already have when you have it).

*Use case.* The same model with other lengths, without rebuilding anything symbolic:

```julia
geometry = RP.Geometry(model)                                              # the URDF lengths
long     = RP.Geometry(model, merge(RP.default_parameters(), (; L2 = 0.18)))   # a longer tibia
```

#### `cartesian_coordinates(q, geometry)`, `swing_foot_position(q, geometry)`

*Description.* For `q = [q₁, q₂, q₃, q₄]` (an `SVector{4}`), the world `(x, z)` positions of
`(grounded_foot, grounded_knee, hip, swing_knee, swing_foot)`, the grounded foot (ankle) at the
origin; `swing_foot_position` is the last of them. Angles are counter-clockwise: the absolute angle
of the grounded tibia is `q₁`, of the grounded femur `q₁ + q₂`, of the swing femur
`φ₃ = q₁ + q₂ + q₃` and of the swing tibia `φ₄ = φ₃ + q₄`; the swing leg hangs downwards.

*Use case.*

```julia
foot = RP.swing_foot_position(SVector(-0.05, 0.0, -0.2, -0.5), geometry)   # [-0.1485, 0.052]
```

### 4.4 The design model

The robot is velocity-controlled: the three motor joints track commanded velocities
`u = [v₂, v₃, v₄]`, held constant over each sampling period (their accelerations are zero), and only
the ankle has non-trivial dynamics. The design model describes the ankle by its **generalized
momentum** `p₁ = ∂L/∂q̇₁ = (D(q) q̇)₁`, in the first slot of the 5-D state, so that the carving, the
target set and the plots are written for `x = [ankle, q₁, q₂, q₃, q₄]`.

The same system can be written with the ankle *velocity*, `x = [q̇₁, q₁…q₄]`, `ẋ = [q̈₁, q̇₁, v₂, v₃, v₄]`
with `q̈₁ = h₁/D₁₁` (the first row of `D q̈ = h + Bτ`, whose row of `B` is zero, with zero actuated
accelerations). That form is not part of the module any more (it served a centre-simulation example
that was removed); `tests/test_robot_problem.jl` keeps it as an *independent* derivation, from the
mass matrix and the bias vector alone, to check the normal-form model against: the RigidBodyDynamics
reference values of the URDF are expressed in it, and the test integrates it and differentiates
`p₁` along the flow to compare with `ṗ₁` of `reduced_dynamics`.

#### `reduced_dynamics(model) -> NamedTuple`

*Description.* Extracts from the normal form the 5-dimensional velocity-controlled model with the
**ankle momentum first**: the state is `x = (p₁, q₁, q₂, q₃, q₄)`, the input `u = v_a = (q̇₂, q̇₃, q̇₄)`,

```
ẋ = (ṗ₁, q̇₁, v₂, v₃, v₄) ,     q̇₁ = D₁₁⁻¹ (p₁ − D₁ₐ v_a) ,     ṗ₁ = ∂L/∂q₁ ,
```

i.e. the `p_u` and `q̇` rows of the normal form. The dynamics of the actuated velocities are dropped:
the commanded velocities are the input, held constant over each sampling period. It checks that no
torque appears. Returns `(; f, x, u)`, all symbolic, with the parameters still free.

*Use case.* A velocity-controlled design model of an underactuated robot, where only the unactuated
momentum has non-trivial dynamics, then the numbers and the local-growth-bound abstraction. This is
what `biped_4d_normal_form.jl` does; on a small domain of 5 cells per axis it gives the abstraction
of 3125 cells, 27 inputs and 157,947 transitions (160,275 with the global growth bound,
`refinements = 0`):

```julia
reduced = RP.reduced_dynamics(model)                  # x = [p₁, q₁…q₄], u = [v₂, v₃, v₄], parameters symbolic
f       = RP.substitute_parameters(model, reduced.f, vals)
disc = RP.default_discretization(; dx = 0.1, dq1 = 0.1, dq̇1 = 0.05, tstep = 0.05)
h    = MP.get_h(disc.state_grid)
X    = LazySets.Hyperrectangle(; low = -2.5 .* h, high = 2.5 .* h)           # 5 cells per axis
U    = LazySets.Hyperrectangle(; low = fill(-disc.u_max, 3), high = fill(disc.u_max, 3))
lgb  = LGA.LocalGrowthBound(f, reduced.x, reduced.u, X, U, disc.tstep)
opt  = LGA.build_abstraction(lgb; state_grid = disc.state_grid, input_grid = disc.input_grid, print_level = 0)
sym  = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
# later: MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem); MOI.optimize!(opt)
```

### 4.5 The footstep scenario

The state is `x = [ankle, q₁, q₂, q₃, q₄]`, where `ankle` is `q̇₁` or `p₁`. Constraints and targets
are Cartesian conditions on the swing foot, pulled back to the grid.

#### `infeasible_cells(X, state_grid, obstacle, geometry; obstacle_name = "obstacle") -> Carving`

*Description.* The cells of the 5-D box `X` to exclude: cells where the swing foot *might* touch the
Cartesian `obstacle` (a `LazySet` or a vector of them; sound, by a Lipschitz bound on the swing-foot
map, see `swing_foot_deviation_bound`) and cells where the swing foot is surely below the ground
(tolerant, so footholds on the ground stay). The decision is taken on the 4-D angle cells and lifted
along the ankle axis. Only the cells lying entirely inside `X` (`MP.INNER`, the default of the
abstraction) are examined.

The constraints are **named**, so that a failure can say which one is responsible: the ground
constraint is `"ground"` (`RP.GROUND_NAME`), and the obstacles are named by `obstacle_name`, one string
for all of them (default `"obstacle"`) or one string per obstacle. The result is a `Carving`: its
field `cells` is the `MP.CellUnion` to remove (what `system` and `UT.set_minus` take), and `causes`
lists `name => carved 4-D angle cells` for every constraint, `"ground"` first. `carved_by(carving, pos)`
gives the names of the constraints that removed a given 5-D cell.

*Use case.* The foot-height constraint and the obstacle in the origin that forces the foot to lift:

```julia
carving = RP.infeasible_cells(X_box, disc.state_grid, [triangle], geometry; obstacle_name = "lift constraint")
removed = carving.cells
X       = UT.set_minus(X_box, removed)
first.(carving.causes)                      # ["ground", "lift constraint"]
```

`swing_foot_gradient_bound(geometry)` and `swing_foot_deviation_bound(geometry, radius)` are the
per-coordinate Lipschitz constants of the swing-foot map and the displacement they allow within a
cell of half-widths `radius`; the latter is the "Lipschitz carving margin" printed by the examples.

#### `target_cells(state_grid, foothold, geometry, X; nsub = 3, velocity = nothing, tolerance = 0.0)`, `target_set(…)`

*Description.* The 5-D cells containing a configuration with the grounded foot at the origin and the
swing foot at `foothold`: `(q₁, q₂)` are swept over their cells (`nsub` samples per axis per cell),
`(q₃, q₄)` follow from the closed-form inverse kinematics of the swing leg (knee bent backwards). The
ankle coordinate is free: every cell of `X` along it, or with `velocity = (vlo, vhi)` every cell lying
entirely inside that range. `target_set` is the `CellUnion` of these cells (it errors if there are
none).

With `tolerance > 0` the foot only has to land *near* the foothold: within `tolerance` metres along `x`
and between the foothold height and `tolerance` above it (never below the ground level of the
foothold), sampled every centimetre (so that the targets of increasing tolerances, multiples of 1 cm,
are nested). **This matters for the synthesis.** With `tolerance = 0`
the target is a thin sliver of the angle space (two constraints, the foot position, on four angles),
and every step the over-approximated ankle spreads over two or three cells (§3, boundary convention):
a controller certified against the worst case cannot land inside a sliver, whatever the start. On the
default grid, 24 combinations of starts, ranges and footholds with `tolerance = 0` all fail, and the
controllable set from a 20 cm foothold only contains starts with steps of about 6 cm; with 3 cm,
starts with a 19 cm step are controllable.

*Use case.* A target region for a reach specification, here with a bounded ankle variable and a 3 cm
tolerance:

```julia
T = RP.target_set(disc.state_grid, SVector(0.2, 0.0), geometry, X_box; velocity = (-0.2, 0.2), tolerance = 0.03)
```

#### `default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1, input_steps = 1)`

*Description.* The 5-D state grid (steps `(dq̇1, dq1, dx, dx, dx)`), the 3-D input grid of step
`du = dx / tstep` and the sampling time, as `(; state_grid, input_grid, tstep, u_max, du, input_steps)` with
`u_max = input_steps · du`. `tstep · du = dx` keeps the three motor coordinates on the lattice (an
input `m · du` moves them by exactly `m` cells); the ankle coordinates follow the nonlinear dynamics
and are rounded to the grid. `dq̇1` is the step of the first coordinate, the ankle momentum `p₁`. Every
motor command is an integer multiple of `du` in `[-input_steps · du, input_steps · du]`:
`2 · input_steps + 1` values per joint, `(2 · input_steps + 1)³` inputs (27 by default); `RP.input_count(disc;
one_joint_per_step = false)` returns that number (`1 + 6 · input_steps` with one moving joint per step).

*Use case.* A finer ankle resolution and a shorter sampling period:

```julia
disc = RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05)
disc2 = RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05, input_steps = 2)  # 125 inputs
```

#### `abstraction_size(X, state_grid, removed = nothing; ninputs = 27, threads = Threads.nthreads())`

*Description.* How big the abstraction will be, **before computing it**: `cells_in_box`, the number of
cells lying entirely inside the box `X` (the cells the abstraction keeps, `MP.INNER`); `cells`, what is
left after removing the carved ones (`removed`: a `Carving` or its `cells`); `pairs = cells · ninputs`,
the number of `(cell, input)` pairs whose image must be computed; and rough estimates
`transitions = RP.TRANSITIONS_PER_PAIR · pairs` (≈ 4 per pair) and `minutes` of wall time on `threads`
threads (`RP.CPU_SECONDS_PER_PAIR`, ≈ 0.2 ms of CPU per pair). The experiment prints it before building
the abstraction.

The cost is linear in `pairs`, and the number of cells is the **product** of the cell counts along the
five axes, so a state box that looks only slightly wider can be much bigger. With the grid of the
experiment (steps `(0.02, 0.05, 0.1, 0.1, 0.1)`): the box `[-0.2, 0.2] × [-0.3, 0.2] × [-0.5, 0.5] ×
[-0.8, 1.0] × [-1.0, 0.2]` has 19 × 9 × 9 × 17 × 11 = 287 793 cells (7.8 M pairs, 34 M transitions, about
3 minutes on 10 threads; 176 225 cells, 4.8 M pairs and 19 M transitions with the obstacle); the box
`[-0.5, 0.5]` on the ankle momentum and on `q₁` and `q₂`, and `[-1, 1]` on `q₃` and `q₄`, has
49 × 19 × 9 × 19 × 19 = 3 024 819 cells, ten times more (1 961 470 after carving the obstacle): 53 M
pairs, about 210 M transitions (several gigabytes of memory) and at least 18 minutes of computation,
much more once memory runs short. From the terminal it looks like the run is stuck after "Number of
forward images", while the abstraction is being computed: there is no progress output.

*Use case.* Check a box before launching:

```julia
sized = RP.abstraction_size(X_box, disc.state_grid, carving)
sized.cells, sized.pairs, round(sized.minutes; digits = 1)      # (287793, 7770411, 3.1) on 10 threads
```

#### `step_problem(concrete_system, state_grid; q0, foothold, geometry, initial_velocity = nothing, target_velocity = nothing, foothold_tolerance = 0.0)`

*Description.* One footstep as a reach specification: from the **configuration**
`q0 = [q₁,q₂,q₃,q₄]` with the first state coordinate (the ankle momentum `p₁`) anywhere in
`initial_velocity`, reach a configuration placing the swing foot at `foothold` (within
`foothold_tolerance` metres, see `target_cells`), with the first coordinate anywhere in
`target_velocity`. Both ranges are `(vlo, vhi)` and default to `nothing`: the whole interval of
definition of the state domain. Nothing is imposed on the ankle at the start or at the end: the swing
foot has no ground reaction force, so requiring rest is needlessly conservative. The initial set is a
box inside the grid cell of `q0` (snapped to the grid) times the initial range: the controller must
work from **every cell of the first coordinate the range touches** (Dionysos lifts the initial set with
the `OUTER` mode: a range `(-0.011, 0.011)` on cells of width `0.02` touches three cells, and
`(-0.031, 0.031)` five). The target is the set of cells lying **entirely** inside the target range
(`INNER` mode). It throws an error if the cell of
`q0` is outside the domain or carved out by the obstacle/ground margin (the initial set would be
empty, which the synthesis can only report as a failure without saying why), or if a range is
reversed, outside the domain, or contains no whole cell.

*Use case.* The footstep that has a controller (no obstacle, see the experiment document): from feet
together, a 20 cm step landing within 3 cm of the foothold, starting with a bounded ankle momentum
and ending with any:

```julia
pb = RP.step_problem(concrete_system, state_grid; q0 = SVector(0.0, -0.3, 0.4, -0.5),
                     foothold = SVector(0.2, 0.0), geometry = geometry,
                     initial_velocity = (-0.031, 0.031),   # target_velocity = nothing: any
                     foothold_tolerance = 0.03)
```

#### `start_state(state_grid, q0, X; velocity = nothing, removed = nothing, geometry)`

*Description.* The concrete start of a closed-loop run: the configuration `q0` with the first state
coordinate (the ankle momentum) in the middle of `velocity` (a `(vlo, vhi)` range clipped to
the box `X`, by default the whole interval of `X`), snapped to the centre of its grid cell so that the
run starts where the abstraction does. With `removed` (the `Carving` of `infeasible_cells`, or just its
`cells`) it errors if the start cell was carved out, *before* any expensive abstraction. With a
`Carving`, the error **names the constraint** that did it and says how the swing foot at the start
violates it:

```
the start cell is carved out by "lift constraint" (the swing foot is within 7.5 cm of it): the swing foot
at the start is at (x, z) = (0.0, 0.0) m; choose another q0 or refine the grid (the margin scales with the
cell sizes)
```

For the ground the reason is "the swing foot is more than 7.5 cm below the ground", and a cell removed
by both constraints names both. (7.5 cm is the Lipschitz carving margin of the grid.) With only the
cells, the message is the generic one without the names.

*Use case.*

```julia
x0 = RP.start_state(disc.state_grid, q0, X_box; velocity = (-0.05, 0.05), removed = carving, geometry)
```

#### `closed_loop(optimizer, problem, x0; horizon = 400, ankle_name = "ankle") -> Trajectory`

*Description.* Simulates the concrete (RK4) system from `x0` under the controller that `optimizer`
synthesized for `problem` (the problem of `step_problem`), for at most `horizon` steps, stopping when
the target set is reached or when the controller has no action at the current state (the concrete
ankle drifted away from the rounded abstract one). It prints the number of steps, whether the target
was reached, whether the run stopped outside the controller's domain, and the ranges of `q₁` and of the
first state coordinate along the run (`ankle_name`, e.g. `"p₁"`). It is the end-to-end check of a
controller and of the abstraction under it. The returned trajectory is what `plot_footstep` draws.

*Use case.* After `MOI.set(optimizer, "concrete_problem", problem)` and `MOI.optimize!(optimizer)`:

```julia
traj = RP.closed_loop(optimizer, problem, x0; ankle_name = "p₁")
```

#### `draw_robot!(fig, x; geometry, alpha)`, `system_plot!(; geometry, obstacle, foothold, xlims, ylims)`

*Description.* `draw_robot!` draws the grounded (blue) and swing (red) legs of a state;
`system_plot!` returns the closure `(fig, x, u) -> fig` that adds the ground, the obstacle and the
foothold, which the trajectory dashboards (`Dionysos.animate_trajectory_dashboard`) call at every
step.

*Use case.* `plot_robot = RP.system_plot!(; geometry, obstacle, foothold)`, to draw a state or to
build your own dashboard.

#### `plot_footstep(traj, state_grid, removed; x_lb, x_ub, tstep, file_prefix, geometry, obstacle, foothold, ankle_label, title, xlims, ylims, snapshot_every = 5, dynamics = nothing, subframes = 1, fps = nothing)`

*Description.* The figures of a closed-loop footstep. From `traj` (a
`Dionysos.System.Trajectory` with inputs, in the state `[ankle, q₁, q₂, q₃, q₄]`) it saves three files
and returns their paths as `(; snapshot, ankle, animation)`:

* `file_prefix * ".png"`: the robot every `snapshot_every` steps, with the ground, the obstacle and
  the foothold;
* `file_prefix * "_ankle.png"`: the ankle angle `q₁` and the ankle coordinate (named by
  `ankle_label`, e.g. `"p₁ (N m s)"`) against time;
* `file_prefix * ".gif"`: the animated dashboard, with the robot, the state panel `(q₃, q₄)` over the
  box `[x_lb, x_ub]` with the slice of the carved cells `removed` at the current `(q₁, q₂)` as a grey
  background (the carving does not depend on the ankle coordinate), and the input panel `(v₃, v₄)`.

It returns `nothing`, drawing nothing, when the trajectory has no input (the controller had no
action at the start).

A closed loop has one state per sampling period, so a run of a few steps is a flash. With `dynamics`
(the continuous dynamics `(x, u) -> ẋ`, e.g. `lgb.dynamics.f`) and `subframes > 1`, each step is
re-integrated (RK4, the input held constant) at `subframes` intermediate instants for the ankle plot and
the animation; `fps` defaults to the value that makes the gif last about four seconds (between 5 and 20).

*Use case.* After a closed-loop run, a smooth animation of a two-step footstep:

```julia
traj = RP.closed_loop(optimizer, problem, x0; ankle_name = "p₁")
RP.plot_footstep(traj, disc.state_grid, removed; x_lb, x_ub, tstep = disc.tstep,
                 geometry, obstacle, foothold, file_prefix = "my_footstep", ankle_label = "p₁ (N m s)",
                 dynamics = lgb.dynamics.f, subframes = 10)
```

### 4.6 What changed in this folder

* **Counter-clockwise convention.** All joint angles, velocities, inputs and the ankle coordinate are
  positive counter-clockwise, like `Multibody2D`: a positive `q₁` leans the robot backwards, a
  positive `q₃` swings the swing foot forwards, and the knee bends backwards for a negative `q₄`.
  (`BipedFourLink.jl` and the URDF keep the clockwise convention of RigidBodyDynamics; see
  `BipedFourLink.md`.)
* **No more RigidBodyDynamics.** The file used to load `RigidBodyDynamics` and `BipedFourLink.jl`
  only to compute the ankle acceleration of the velocity-controlled model. `BipedFourLink.jl` is no
  longer part of the design chain: it only contains simulation and visualization.
* **Footstep problem with any ankle momentum.** `balanced_start` (an equilibrium at rest) was removed;
  the initial and target ranges of the ankle coordinate are free (`step_problem`, above).
* **Consolidation: nothing is defined twice.** `BipedModel.jl` and `UnderactuatedLocalAbstraction.jl`
  described the same four-link model and the same 5-D systems as `robot_problem.jl`, with a second
  copy of the lengths, of the kinematics and of the discretization. They were merged into it.
* **Only the normal-form experiment is kept.** The centre-simulation example
  (`biped_4d_underactuated.jl`, the ankle described by its velocity, not a certified abstraction) was
  removed with what only it used: `VelocityControlledDynamics`, `RobotDomainConfig` and `system`. The
  velocity form of the dynamics survives in the tests as an independent check (§4.4).
* **Pre-flight size, explained start errors, foothold tolerance.** `abstraction_size` prints the cost of
  a box before the abstraction starts; the carving constraints are named (`Carving`) so that the start
  error says which one removed the start cell; `target_cells` and `step_problem` accept a foothold
  tolerance, without which no controller is found (§4.5).

| Removed | Where it went |
| :--- | :--- |
| `BipedModel.jl`: `symbolic_model`, `default_parameters` | `robot_problem.jl`, unchanged |
| `BipedModel.ankle_acceleration(; values)`, `VelocityControlledDynamics`, `RobotDomainConfig`, `system` | removed with the centre-simulation example; the velocity form of `q̈₁ = h₁/D₁₁` is a helper of `tests/test_robot_problem.jl` |
| `Geometry(; tibia, femur)`, the closed-form `cartesian_coordinates` and `check_kinematics` (a check of one copy against the other) | `Geometry(model, values)`, **compiled from the symbolic positions**; the hand-derived closed form lives only in the tests, as an independent check |
| `UnderactuatedLocalAbstraction.numerical_dynamics(model, reduced, values)` | `substitute_parameters(model, expressions, values)`, shared with the kinematics |
| `UnderactuatedLocalAbstraction.reduced_dynamics` (state `(q, p₁)`) and the permutation `[5, 1, 2, 3, 4]` in the experiment | `reduced_dynamics`, with `p₁` first |
| `default_discretization`, `local_growth_bound`, `build_local_abstraction` of `UnderactuatedLocalAbstraction` | `RP.default_discretization` and `LGA.LocalGrowthBound` / `LGA.build_abstraction`, which the experiment calls directly (the snippet of §4.4 is the whole chain) |
| `check_reduced_dynamics`, `check_abstraction`, `main` of `UnderactuatedLocalAbstraction` | tests (`tests/test_robot_problem.jl`): equivalence of the velocity and momentum forms, soundness by sampling and local ⊆ global on a small grid |

`BipedFourLink.jl` is **not** merged: it is the only file that needs RigidBodyDynamics and MeshCat
(heavy dependencies useless for design), it reads the URDF that MeshCat draws, takes motor torques
as inputs and uses the clockwise convention. Its consistency with the design model is enforced by
`tests/test_biped_four_link.jl` instead. Inside it, the tibia length of the foot point (three copies)
became one table, and the two state setters (locked ankle, compliant feet) one function.
`Multibody2D.jl` and `LocalGrowthAbstraction.jl` are kept apart because they are generic and are
tested on other systems (double pendulum, cart-pole, acrobot, damped pendulum).

---

## 5. Tests and known limitations

Tests (run from the repository root, with an environment that has the packages of §1.2):

```
julia --project=<env> problems/BipedRobot/4D_model_underactuated/tests/test_multibody2d.jl                 # 172 checks
julia --project=<env> problems/BipedRobot/4D_model_underactuated/tests/test_local_growth_abstraction.jl    # 5573 checks
julia --project=<env> problems/BipedRobot/4D_model_underactuated/tests/test_robot_problem.jl               # 200 checks
julia --project=<env> problems/BipedRobot/4D_model_underactuated/tests/test_biped_four_link.jl             # 53 checks; needs RigidBodyDynamics and MeshCat
```

`test_multibody2d.jl` checks the energies and the normal form against closed-form double-pendulum
and cart-pole results and, numerically, against the second-order dynamics `D q̈ = h + Bτ`
(finite differences for `Ḋ`; one, two and three unactuated coordinates, and fully actuated).
`test_local_growth_abstraction.jl` checks the Jacobian bound against Dionysos'
`compute_jacobian_bound`, the flow enclosure against simulated tubes, soundness by sampling,
`local ⊆ global` transitions, the fallback, the boundary tolerance, and a full run on the acrobot
normal form from `Multibody2D`. `test_robot_problem.jl` checks `robot_problem.jl`:

* the dynamics: the velocity form of `q̈₁` (a helper of the test file) against reference values
  computed once with RigidBodyDynamics on the URDF, a static equilibrium, and the equivalence of the
  velocity and momentum models (integrating one and differentiating `p₁` along the flow);
* the kinematics against a hand-derived closed form (for other lengths too);
* the local-growth-bound abstraction of the normal-form model on a 3⁵-cell grid (soundness by
  sampling, local ⊆ global), and on the same abstraction `start_state` and `closed_loop`: a controller
  is synthesized and the closed loop reaches the target;
* the scenario: sound carving and its named constraints, the start-state error that names the
  constraint (obstacle, ground, both), inverse kinematics, target cells and the foothold tolerance,
  `step_problem` and its errors, the size of the abstraction (`abstraction_size`, against a brute-force
  count and for the boxes of the experiment), and `plot_footstep` on a synthetic trajectory, with and
  without re-integration of the steps.

`test_biped_four_link.jl` checks the simulation (state round trips for both mechanisms, energy
conservation, compliant contact, reaction force, plots, MeshCat) and its agreement with the design model
(torque-level accelerations, forward kinematics).

Known limitations:

* the endpoint arithmetic of the flow enclosure is not outward rounded;
* the local growth bound needs a short enough sampling period (§3.5);
* the gain over the global growth bound is modest when the Jacobian varies little over the domain
  (1.5 % on the 5⁵ biped grid of §4.4, 7 % on a 9⁵ grid, 34 % on the acrobot test). The image of
  every enabled `(cell, input)` pair is inside its global image, so the local abstraction enables at
  least the pairs of the global one with fewer successors; on a tiny domain it can still have more
  transitions in total, because pairs whose global image leaves the domain become enabled (3⁵ biped
  grid: 4151 against 4135).
