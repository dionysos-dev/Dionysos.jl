# `BipedFourLink`: simulation and visualization of the four-link biped

`BipedFourLink.jl` simulates and visualizes the planar four-link biped with
[RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl), Plots and MeshCat,
from `four_link_biped.urdf`. It is the tool to use once a controller exists and the **whole
robot** is to be simulated with **motor torques** as inputs.

It is *not* the design model. The mass matrix, bias terms, velocity-controlled model and normal form
used to design the controller are in `Multibody2D.jl` / `robot_problem.jl` (symbolic, no
RigidBodyDynamics); see `symbolic_multibody_and_local_abstraction.md`.

```julia
include("BipedFourLink.jl")
using .BipedFourLink

model = build_model("four_link_biped.urdf")      # a SimulationModel
```

Dependencies: `RigidBodyDynamics`, `StaticArrays`, `Plots`, `MeshCatMechanisms` (and `MeshCat`).

## Contents

1. [Coordinates and the angle convention](#1-coordinates-and-the-angle-convention)
2. [Relation to the design model](#2-relation-to-the-design-model)
3. [The two mechanisms: locked ankle and compliant feet](#3-the-two-mechanisms-locked-ankle-and-compliant-feet)
4. [Simulating with torques](#4-simulating-with-torques)
5. [Analysing a simulation: reactions, friction, slip](#5-analysing-a-simulation-reactions-friction-slip)
6. [Plots and MeshCat](#6-plots-and-meshcat)
7. [Reference](#7-reference)

---

## 1. Coordinates and the angle convention

| Coordinate | Meaning | Actuation |
| --- | --- | --- |
| `q₁` | grounded ankle relative to the world vertical | **none** (unactuated) |
| `q₂` | grounded femur relative to grounded tibia (grounded knee) | direct torque `τ₂` |
| `q₃` | swing femur relative to grounded femur (hip) | direct torque `τ₃` |
| `q₄` | swing tibia relative to swing femur (swing knee) | direct torque `τ₄` |

```julia
q  = [q1, q2, q3, q4]
q̇ = [q̇1, q̇2, q̇3, q̇4]
u  = [τ2, τ3, τ4]            # the inputs: torques on the three motor joints
```

* `q₁ = 0`: grounded tibia vertical; `q₂ = 0`: grounded femur aligned with the tibia; `q₃ = 0`: the
  two femurs side by side; `q₄ = 0`: the swing leg straight.
* The grounded tibia and femur are drawn from their proximal joint along local `+z`; the swing
  femur and tibia along local `-z`. There is no torso.

> **Angle convention: clockwise.** Angles, velocities and torques follow the URDF: the joints rotate
> about `+y`, so a positive `q₁` leans the robot towards `+x`. **`Multibody2D` and
> `robot_problem.jl` are counter-clockwise.** Going from one to the other, `q`, `q̇`, `q̈` and `τ`
> change sign; world positions do not.

## 2. Relation to the design model

The design model (`UnderactuatedRobotProblem.symbolic_model`, counter-clockwise; `RP` below) and
this simulation describe the **same robot**. `four_link_biped.urdf` has the same masses and link
lengths as `RP.default_parameters()`, uniform
thin rods (`I = m L²/12` about their centre of mass, the `<inertial>` origin at the middle of each
link) and 0.4 kg point-mass motors with no rotor inertia (two at the hip). Consequently their
torque-level dynamics agree: for any state and torque,

```
accelerations(model, q, q̇, u)  ==  -( D(−q) \ ( h(−q, −q̇) + B·(−u) ) ) ,
```

with `D`, `h` the mass matrix and bias vector of `Multibody2D` and `B = [0 0 0; I₃]` (the signs are
those of the clockwise/counter-clockwise conversion of §1). This is checked in
`tests/test_biped_four_link.jl`, together with the forward kinematics (joint and foot positions of
the URDF against `RP.Geometry`, with `q_design = -q_urdf`), and the ankle acceleration of the
velocity-controlled model against reference values of this simulation in
`tests/test_robot_problem.jl`. A controller designed on the symbolic model is therefore simulated
here on the same inertias and geometry.

## 3. The two mechanisms: locked ankle and compliant feet

**Locked ankle** (`model.mechanism`, the default). `q₁` is a revolute joint from the world to the
grounded tibia: the grounded foot is pinned at the origin, an exact, perfectly rigid ankle support.
It is a four-DoF tree, so the ankle reaction can be reconstructed from inverse dynamics without
floating-base coordinates or loop multipliers.

```julia
state = state_from_coordinates(model, q, q̇)
q_again, q̇_again = coordinates_from_state(model, state)

p_knee       = world_position(model, state, "grounded_knee")
p_swing_foot = foot_position(model, state; foot = :swing)      # [x, y, z]; z < 0 is below the ground
```

**Compliant feet** (built on demand). A six-DoF planar tree with coordinates
`[ankle x, ankle z, q₁, q₂, q₃, q₄]`. It uses a planar root joint, so it cannot move laterally, roll
or yaw. The grounded ankle point and the swing-tibia endpoint are both subject to a Hunt–Crossley
normal contact and Coulomb-limited friction (`μ`; `normal_stiffness`, `normal_damping`,
`tangential_stiffness`, `tangential_damping`). Raising a stiffness requires reducing `Δt`.

## 4. Simulating with torques

`controller!(u, t, q, q̇)` fills the three motor torques `u = [τ₂, τ₃, τ₄]` at every integration
step from the current coordinates (`[q₁…q₄]`, `[q̇₁…q̇₄]`):

```julia
function controller!(u, t, q, q̇)
    u .= 0.0                                   # free fall of the unactuated robot
end

# exact anchored ankle
ts, qhist, q̇hist, qs, vs = simulate_trajectory(model, q0, q̇0, 2.0, controller!;
                                               foot_model = :locked, Δt = 1e-3)

# compliant two-foot contact
ts, qhist, q̇hist, qs, vs, status = simulate_trajectory(model, q0, q̇0, 2.0, controller!;
                                               foot_model = :compliant, μ = 0.8, Δt = 1e-4,
                                               return_contact_status = true, contact_foot = :both)
```

`qhist`, `q̇hist` always contain only `[q₁, q₂, q₃, q₄]`; the full RigidBodyDynamics configurations
`qs` and velocities `vs` (with the ankle translation for the compliant model) are kept for animation
and contact analysis.

**An inner velocity loop** that tracks reference actuated velocities coming from the design model is
written in the same hook, converting signs at the boundary (the references come from a
counter-clockwise model):

```julia
v_ref_ccw = [0.5, -0.5, 0.5]                     # [v₂, v₃, v₄] from the design model
function velocity_loop!(u, t, q, q̇)
    τ_ccw = Kp .* (v_ref_ccw .- (-q̇[2:4]))       # q̇ is clockwise here
    u .= -τ_ccw                                  # the simulator expects clockwise torques
end
```

A caution found while testing this: the actuated joints have *very small* effective inertias (down to
`1.7e-4 kg m²` for the swing knee), and `simulate_trajectory` integrates explicitly. A proportional
gain with `Kp / I ≳ 1/Δt` diverges (`Kp = 2` with `Δt = 1e-3` gave `NaN`); choose the gain and `Δt`
together (a `Kp = 0.01` loop is stable at `Δt = 1e-3` but far too weak to reject gravity and the
falling ankle, so a real controller needs a stiffer loop with a smaller step or an implicit/
feed-forward design).

## 5. Analysing a simulation: reactions, friction, slip

For the locked model, the world reaction at the grounded ankle, and a Coulomb friction-cone test of
it (the anchored-ankle assumption is physically feasible only if the ground can supply it):

```julia
q̈     = accelerations(model, q, q̇, u)                 # forward dynamics under torques u = [τ₂, τ₃, τ₄]
force = ankle_reaction_force(model, q, q̇, q̈)         # [Fx, Fy, Fz], world frame
wrench = ankle_reaction_wrench(model, q, q̇, q̈)       # the full wrench
ok    = in_friction_cone(force, 0.8; normal_axis = 3, normal_sign = 1)    # ‖F_t‖ ≤ μ F_n and F_n ≥ 0
```

For a compliant trajectory, `status.grounded` and `status.swing` each report
`(was_in_contact, stayed_in_place, slipped, max_slip, contact_samples)`; the slip test uses the
tangential world-`x` displacement during each continuous contact phase (`slip_tolerance`, default
`1e-4` m). A saved trajectory can be analysed separately:

```julia
status = contact_slip_status(model, qs; contact_foot = :grounded)
```

## 6. Plots and MeshCat

```julia
forces = [ankle_reaction_force(model, qhist[k], q̇hist[k],
              accelerations(model, qhist[k], q̇hist[k], u_of(ts[k]))) for k in eachindex(ts)]
plots = plot_trajectory(ts, qhist, q̇hist; forces)     # (positions, velocities, forces) Plots figures
savefig(plots.positions, "positions.png")

vis = animate_trajectory(model, ts, qs; foot_model = :locked)       # MeshCat; open(vis) in a browser
```

For a compliant trajectory pass `foot_model = :compliant` to `animate_trajectory`. Soft-contact forces
are internal contact-model forces; use the locked model's `ankle_reaction_force` for pointwise
friction-cone feasibility tests.

## 7. Reference

| Name | Description |
| :--- | :--- |
| `SimulationModel`, `build_model(urdf_file)` | the locked-ankle RigidBodyDynamics mechanism built from a URDF (the struct was formerly called `BipedModel`) |
| `state_from_coordinates(!)`, `coordinates_from_state` | `[q, q̇]` ↔ RigidBodyDynamics state, for either mechanism (`state_from_coordinates(model, q, q̇; mechanism)` takes the mechanism of the state to build, the locked one by default) |
| `world_position`, `foot_position` | world position of a URDF joint/body, of the grounded/swing foot |
| `accelerations(model, q, q̇, u)` | forward dynamics under motor torques `u = [τ₂, τ₃, τ₄]` (locked model) |
| `ankle_reaction_wrench`, `ankle_reaction_force`, `in_friction_cone` | anchored-ankle reaction and friction-cone test |
| `simulate_trajectory(model, q0, q̇0, T, controller!; foot_model, μ, Δt, …)` | locked or compliant simulation with a torque controller |
| `contact_slip_status(model, qs; …)` | slip analysis of a compliant trajectory |
| `plot_trajectory`, `animate_trajectory` | Plots figures; MeshCat animation |
| `COORDINATE_NAMES`, `INPUT_NAMES` | labels |

**Removed** (the design-model dynamics, now in `Multibody2D` / `robot_problem.jl`): `coordinate_mass_matrix`,
`bias_forces`, `coordinate_inverse_dynamics`, `velocity_controlled_derivative`,
`velocity_controlled_contact_estimate`. Their replacements: `Multibody2D.mass_matrix`,
`Multibody2D.bias_vector`, `normal_form`, `RP.reduced_dynamics`.

Tests: `julia --project=<env> tests/test_biped_four_link.jl` (53 checks: state round trip for both
mechanisms, forward kinematics and torque-level agreement with the design model, energy
conservation, compliant contact, reaction force and friction cone, plots, MeshCat visualizer).
