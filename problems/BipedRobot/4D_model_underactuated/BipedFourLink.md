# `BipedFourLink`: four-coordinate planar biped

`BipedFourLink.jl` and `four_link_biped.urdf` are the replacement package for
the torso-free biped. Copy both files into the toolbox source directory.

```julia
include("BipedFourLink.jl")
using .BipedFourLink

model = build_model("four_link_biped.urdf")
```

Dependencies are `RigidBodyDynamics`, `StaticArrays`, `Plots`, and
`MeshCatMechanisms`.

## Coordinate convention

The locked design model has exactly four positions and four velocities:

| Coordinate | Meaning | Actuation |
| --- | --- | --- |
| `q₁` | grounded ankle relative to the world vertical | underactuated |
| `q₂` | grounded femur relative to grounded tibia | direct torque `τ₂` |
| `q₃` | swing femur relative to grounded femur | direct torque `τ₃` |
| `q₄` | swing tibia relative to swing femur | direct torque `τ₄` |

```julia
q  = [q1, q2, q3, q4]
q̇ = [q̇1, q̇2, q̇3, q̇4]
u  = [τ2, τ3, τ4]
```

The grounded tibia and femur are drawn from their proximal joint along local
`+z`. The swing femur and tibia are drawn along local `-z`.

- `q₁ = 0`: grounded tibia is vertical.
- `q₂ = 0`: grounded femur and grounded tibia are aligned.
- `q₃ = 0`: the two femurs are side-by-side.
- `q₄ = 0`: the swing leg is straight.

There is no torso. Two 0.4 kg hip-motor bodies are located at the hip, and a
0.4 kg knee-motor body is located at each knee. They retain the previous motor
inertias but are fixed bodies: the actuator torques act directly on `q₂`,
`q₃`, and `q₄`.

## Locked-ankle design model

The default model roots `q₁` in a world-to-grounded-tibia revolute joint. It
is a four-DoF tree, so the ankle reaction can be reconstructed directly from
inverse dynamics without floating-base coordinates or loop multipliers.

```julia
state = state_from_coordinates(model, q, q̇)
q_again, q̇_again = coordinates_from_state(model, state)

p_knee = world_position(model, state, "grounded_knee")
p_swing_foot = foot_position(model, state; foot = :swing)
```

`foot_position` accepts `foot = :grounded` or `:swing`. Use their third
component to reject planning configurations that go below the ground plane.

## Dynamics and inverse dynamics

All quantities below are already in the four requested coordinates:

```julia
M = coordinate_mass_matrix(model, q, q̇)  # 4 × 4
c = bias_forces(model, q, q̇)             # gravity and velocity terms
q̈ = accelerations(model, q, q̇, u)         # u = [τ₂, τ₃, τ₄]

required = coordinate_inverse_dynamics(model, q, q̇, q̈)
# required.u            -> [τ₂, τ₃, τ₄]
# required.ankle_torque -> required external q₁ torque
```

The physical equation is

```math
M(q)\ddot q + c(q,\dot q) = [\tau_1,\tau_2,\tau_3,\tau_4]^T.
```

For the unactuated design model, `τ₁ = 0`. Therefore a prescribed
acceleration is dynamically feasible only when
`required.ankle_torque ≈ 0`.

## Velocity-controlled design dynamics

For graph-based controller design, the module also supports the kinematic
actuation abstraction in which the three motor-joint velocities are direct
inputs. Its state and input are

```julia
x = [q̇1, q1, q2, q3, q4]
v = [v2, v3, v4]
```

Use:

```julia
ẋ = velocity_controlled_derivative(model, x, v)
```

The returned dynamics are

```math
\begin{aligned}
\dot q_1 &= x_1,\\
\ddot q_1 &= -c_1(q,[x_1,v]^T)/M_{11}(q),\\
\dot q_2 &= v_2,\\
\dot q_3 &= v_3,\\
\dot q_4 &= v_4.
\end{aligned}
```

This is the first, unactuated Euler--Lagrange equation with the direct-velocity
commands treated as piecewise constant over a graph transition. If a planner
also specifies their derivative, use:

```julia
ẋ = velocity_controlled_derivative(model, x, v;
    input_accelerations = [v̇2, v̇3, v̇4])
```

The exact more general formula is

```math
\ddot q_1 = -\frac{M_{1a}\dot v + c_1}{M_{11}}.
```

This distinction matters: instantaneous velocities alone do not determine the
exact unactuated acceleration unless the velocity command is held constant or
its derivative is supplied.

To estimate the physical motor demands and grounded-ankle reaction after a
velocity-controlled transition, use:

```julia
estimate = velocity_controlled_contact_estimate(model, x, v)
estimate.force       # [Fx, Fy, Fz]
estimate.wrench      # full ankle wrench
estimate.q̈           # inferred [q̈1, v̇2, v̇3, v̇4]
estimate.required.u  # physical [τ2, τ3, τ4]
```

This estimate is deliberately a posteriori: the graph dynamics remain the
five-state, three-velocity-input model above.

## Anchored-ankle reaction and friction cone

```julia
wrench = ankle_reaction_wrench(model, q, q̇, q̈)
force = ankle_reaction_force(model, q, q̇, q̈)  # [Fx, Fy, Fz], world frame

feasible = in_friction_cone(force, 0.8; normal_axis = 3, normal_sign = 1)
```

The wrench is the world reaction at the grounded ankle. The friction test
checks `||F_t|| ≤ μ F_n` and `F_n ≥ 0`. Its configurable `normal_sign` makes
the sign convention explicit for integration with other toolbox components.

## Simulating locked or compliant feet

Provide a controller that fills the direct three-motor input:

```julia
function controller!(u, t, q, q̇)
    u .= 0.0
end
```

### Exact anchored ankle

```julia
ts, qhist, q̇hist, qs, vs = simulate_trajectory(
    model, q, q̇, 2.0, controller!;
    foot_model = :locked,
    Δt = 1e-3)
```

### Compliant two-foot contact

```julia
ts, qhist, q̇hist, qs, vs, status = simulate_trajectory(
    model, q, q̇, 2.0, controller!;
    foot_model = :compliant,
    μ = 0.8,
    Δt = 1e-4,
    return_contact_status = true,
    contact_foot = :both)
```

The compliant mechanism is a **six-DoF planar tree**:

```text
[grounded-ankle x, grounded-ankle z, q₁, q₂, q₃, q₄]
```

It uses a planar root joint, so it cannot move laterally, roll, or yaw. The
grounded ankle point and the swing-tibia endpoint are both subject to a
Hunt--Crossley normal contact and Coulomb-limited friction. `μ` controls the
Coulomb coefficient. `normal_stiffness`, `normal_damping`,
`tangential_stiffness`, and `tangential_damping` are also available keyword
parameters. Raising stiffness requires reducing `Δt`.

`qhist` and `q̇hist` always contain only `[q₁, q₂, q₃, q₄]`; the additional
grounded-ankle translation is retained in full `qs` and `vs` for animation or
further contact analysis.

`status.grounded` and `status.swing` each report:

```julia
(was_in_contact, stayed_in_place, slipped, max_slip, contact_samples)
```

The slip test uses tangential world-`x` displacement during each continuous
contact phase. Tune `slip_tolerance` (default `1e-4` m) to match the desired
stiction accuracy. To analyse a saved compliant trajectory separately:

```julia
status = contact_slip_status(model, qs; contact_foot = :grounded)
```

## Plots and MeshCat

For the locked trajectory, reconstruct and plot ankle force:

```julia
forces = Vector{Vector{Float64}}(undef, length(ts))
for k in eachindex(ts)
    u = zeros(3)
    controller!(u, ts[k], qhist[k], q̇hist[k])
    q̈ = accelerations(model, qhist[k], q̇hist[k], u)
    forces[k] = ankle_reaction_force(model, qhist[k], q̇hist[k], q̈)
end

plots = plot_trajectory(ts, qhist, q̇hist; forces)
savefig(plots.positions, "positions.png")
savefig(plots.velocities, "velocities.png")
savefig(plots.forces, "ankle_force.png")

vis = animate_trajectory(model, ts, qs; foot_model = :locked)
```

For a compliant trajectory, pass `foot_model = :compliant` to
`animate_trajectory`. Soft-contact forces are internal contact-model forces;
use the locked model’s `ankle_reaction_force` for pointwise design and
friction-cone feasibility tests.
