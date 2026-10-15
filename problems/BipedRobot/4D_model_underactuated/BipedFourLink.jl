"""
    BipedFourLink

**Simulation and visualization** of the planar four-link biped with RigidBodyDynamics.jl,
MeshCat and Plots, from `four_link_biped.urdf`. The design dynamics (mass matrix, bias terms,
velocity-controlled model, normal form) live in `Multibody2D.jl` / `robot_problem.jl`; this module is
what you use once a controller exists and the whole robot is to be simulated with **motor torques**
as inputs.

Two mechanisms are built from the URDF:

* the **locked-ankle** four-DoF model (`model.mechanism`): the grounded foot is pinned to the origin
  by a revolute joint `q₁` (an exact, perfectly rigid ankle support);
* a **compliant-foot** six-DoF planar model (built on demand): the ankle can translate, and both feet
  are subject to Hunt–Crossley normal contact and Coulomb friction on the ground.

!!! note "Angle convention"
    Angles, velocities and torques here follow the URDF: **positive clockwise** in the `(x, z)`
    plane (a rotation about `+y`; positive `q₁` leans the robot towards `+x`). `Multibody2D` and
    `robot_problem.jl` are **counter-clockwise**, so `q`, `q̇` and `τ` must change sign when passing
    between them (world positions do not).
"""
module BipedFourLink

using StaticArrays
using RigidBodyDynamics
using Plots
using MeshCatMechanisms

const RBD = RigidBodyDynamics

export SimulationModel,
    build_model,
    state_from_coordinates!,
    state_from_coordinates,
    coordinates_from_state,
    world_position,
    foot_position,
    accelerations,
    ankle_reaction_wrench,
    ankle_reaction_force,
    in_friction_cone,
    simulate_trajectory,
    contact_slip_status,
    plot_trajectory,
    animate_trajectory,
    COORDINATE_NAMES,
    INPUT_NAMES

"""
    SimulationModel

The RigidBodyDynamics description of the biped built by [`build_model`](@ref).

`mechanism` is the locked-ankle model (an exact revolute ankle support); `q1_joint` is its tree
joint. The compliant-foot model is constructed separately, with a planar root, when needed.
"""
struct SimulationModel
    mechanism::Any
    q1_joint::Any
    joint_names::NTuple{3, String}
    urdf_file::String
end

const COORDINATE_NAMES = ("q1_ankle", "q2_grounded_knee", "q3_hip", "q4_swing_knee")
const INPUT_NAMES = ("τ2_grounded_knee", "τ3_hip", "τ4_swing_knee")

# The foot of each leg: the distal end of its tibia, as (body, z in the body frame). The grounded
# tibia's frame is at the ankle; the swing tibia's foot is one tibia length (`L2 = 0.172` m in
# `four_link_biped.urdf`) below its frame.
const FEET = (grounded = ("grounded_tibia", 0.0), swing = ("swing_tibia", -0.172))

function _foot(mechanism, foot::Symbol)
    foot in keys(FEET) || throw(ArgumentError("foot must be :grounded or :swing"))
    body_name, z = FEET[foot]
    body = RBD.findbody(mechanism, body_name)
    return body, RBD.Point3D(RBD.default_frame(body), 0.0, 0.0, z)
end

_parse(urdf_file, root_joint_type) = RBD.URDF.parse_urdf(
    urdf_file;
    root_joint_type,
    remove_fixed_tree_joints = false,
    gravity = SVector(0.0, 0.0, -9.81),
)

"""
    build_model(urdf_file)

Create the four-DoF locked-ankle mechanism. `q1` is a ground-ankle revolute joint;
`q2`, `q3`, and `q4` are the grounded knee, hip, and swing knee.
"""
function build_model(urdf_file::AbstractString)
    mechanism = _parse(urdf_file, RBD.Revolute(SVector(0.0, 1.0, 0.0)))
    return SimulationModel(
        mechanism,
        RBD.findjoint(mechanism, "grounded_tibia_to_world"),
        ("grounded_knee", "hip", "swing_knee"),
        String(urdf_file),
    )
end

_joint(model, index, mechanism = model.mechanism) =
    RBD.findjoint(mechanism, model.joint_names[index])

"""
    state_from_coordinates!(state, model, q, q̇ = zeros(4))

Set `state`, a state of either mechanism of `model` (locked ankle, or compliant feet whose ankle
starts at the origin), from `q = [q1, q2, q3, q4]` and `q̇`.
"""
function state_from_coordinates!(
    state,
    model::SimulationModel,
    q::AbstractVector,
    q̇::AbstractVector = zeros(4),
)
    length(q) == 4 || throw(DimensionMismatch("q must have four coordinates"))
    length(q̇) == 4 || throw(DimensionMismatch("q̇ must have four coordinates"))
    mechanism = state.mechanism
    if mechanism === model.mechanism
        RBD.set_configuration!(state, model.q1_joint, q[1])
        RBD.set_velocity!(state, model.q1_joint, q̇[1])
    else
        # Planar root: base coordinates [x, -z, q1].
        root = RBD.findjoint(mechanism, "grounded_tibia_to_world")
        RBD.set_configuration!(state, root, [0.0, 0.0, q[1]])
        RBD.set_velocity!(state, root, [0.0, 0.0, q̇[1]])
    end
    for i in 1:3
        RBD.set_configuration!(state, _joint(model, i, mechanism), q[i + 1])
        RBD.set_velocity!(state, _joint(model, i, mechanism), q̇[i + 1])
    end
    return state
end

"""
    state_from_coordinates(model, q, q̇ = zeros(4); mechanism = model.mechanism)

A new state of `mechanism` (the locked-ankle one by default) set from `q` and `q̇`.
"""
function state_from_coordinates(
    model::SimulationModel,
    q::AbstractVector,
    q̇::AbstractVector = zeros(4);
    mechanism = model.mechanism,
)
    return state_from_coordinates!(RBD.MechanismState(mechanism), model, q, q̇)
end

"""Extract `[q1, q2, q3, q4]` and its derivative from either internal model."""
function coordinates_from_state(model::SimulationModel, state)
    mechanism = state.mechanism
    q = RBD.configuration(state)
    v = RBD.velocity(state)
    if mechanism === model.mechanism
        q1_joint = model.q1_joint
        q1 = q[RBD.configuration_range(state, q1_joint)][1]
        q̇1 = v[RBD.velocity_range(state, q1_joint)][1]
    else
        root = RBD.findjoint(mechanism, "grounded_tibia_to_world")
        qr = RBD.configuration_range(state, root)
        vr = RBD.velocity_range(state, root)
        q1 = q[qr][3]      # Planar root: [x, -z, q1]
        q̇1 = v[vr][3]
    end
    qs = [
        q1;
        [q[RBD.configuration_range(state, _joint(model, i, mechanism))][1] for i in 1:3]
    ]
    q̇s = [q̇1; [v[RBD.velocity_range(state, _joint(model, i, mechanism))[1]] for i in 1:3]]
    return qs, q̇s
end

"""World location of a URDF joint or body."""
function world_position(model::SimulationModel, state, name::AbstractString)
    mechanism = state.mechanism
    frame = try
        RBD.frame_after(RBD.findjoint(mechanism, name))
    catch
        RBD.default_frame(RBD.findbody(mechanism, name))
    end
    return RBD.translation(RBD.transform_to_root(state, frame))
end

"""World position of the grounded or swing foot point (`foot = :grounded` or `:swing`)."""
function foot_position(model::SimulationModel, state; foot::Symbol = :grounded)
    mechanism = state.mechanism
    _, point = _foot(mechanism, foot)
    return RBD.transform(state, point, RBD.root_frame(mechanism)).v
end

function _full_torques(state, model, u)
    length(u) == 3 || throw(DimensionMismatch("u must be [τ2, τ3, τ4]"))
    τ = zeros(eltype(RBD.configuration(state)), RBD.num_velocities(state.mechanism))
    return _set_full_torques!(τ, state, model, u)
end

function _set_full_torques!(τ, state, model, u)
    length(u) == 3 || throw(DimensionMismatch("u must be [τ2, τ3, τ4]"))
    fill!(τ, zero(eltype(τ)))
    for i in 1:3
        τ[RBD.velocity_range(state, _joint(model, i, state.mechanism))[1]] = u[i]
    end
    return τ
end

"""Compute `q̈` under direct motor torques `u = [τ2, τ3, τ4]`; q1 is unactuated."""
function accelerations(model::SimulationModel, q, q̇, u)
    state = state_from_coordinates(model, q, q̇)
    result = RBD.DynamicsResult(model.mechanism)
    RBD.dynamics!(result, state, _full_torques(state, model, u))
    return collect(result.v̇)
end

"""World-frame wrench applied by the world at the grounded ankle for `(q, q̇, q̈)`."""
function ankle_reaction_wrench(model::SimulationModel, q, q̇, q̈)
    state = state_from_coordinates(model, q, q̇)
    result = RBD.DynamicsResult(model.mechanism)
    copyto!(result.v̇, q̈)
    RBD.inverse_dynamics!(
        result.dynamicsbias,
        result.jointwrenches,
        result.accelerations,
        state,
        result.v̇,
    )
    grounded_tibia = RBD.findbody(model.mechanism, "grounded_tibia")
    wrench = RBD.joint_wrench(result, grounded_tibia)
    return RBD.transform(state, wrench, RBD.root_frame(model.mechanism))
end

"""World-frame `[Fx, Fy, Fz]` grounded-ankle reaction force."""
ankle_reaction_force(model::SimulationModel, q, q̇, q̈) =
    collect(RBD.linear(ankle_reaction_wrench(model, q, q̇, q̈)))

"""Test `||F_t|| ≤ μ F_n` and `F_n ≥ 0` for a world-frame contact force."""
function in_friction_cone(force, μ; normal_axis = 3, normal_sign = 1)
    fn = normal_sign * force[normal_axis]
    ft2 = sum(force[i]^2 for i in eachindex(force) if i != normal_axis)
    return fn >= 0 && ft2 <= (μ * fn)^2
end

function _compliant_mechanism(
    model::SimulationModel,
    μ;
    normal_stiffness = 1e5,
    normal_damping = 1.0,
    tangential_stiffness = 5e4,
    tangential_damping = 1e4,
)
    # +x and -z translations give a +y pitch axis: x × (-z) = +y.
    mechanism =
        _parse(model.urdf_file, RBD.Planar(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, -1.0)))
    world = RBD.root_frame(mechanism)
    RBD.add_environment_primitive!(
        mechanism,
        RBD.HalfSpace3D(
            RBD.Point3D(world, 0.0, 0.0, 0.0),
            RBD.FreeVector3D(world, 0.0, 0.0, 1.0),
        ),
    )
    normal = RBD.hunt_crossley_hertz(; k = normal_stiffness, α = normal_damping)
    friction = RBD.ViscoelasticCoulombModel(μ, tangential_stiffness, tangential_damping)
    contact = RBD.SoftContactModel(normal, friction)
    for foot in keys(FEET)
        body, point = _foot(mechanism, foot)
        RBD.add_contact_point!(body, RBD.ContactPoint(point, contact))
    end
    return mechanism
end

function _single_foot_slip_status(
    mechanism,
    qs,
    foot;
    ground_tolerance = 1e-5,
    slip_tolerance = 1e-4,
)
    state = RBD.MechanismState(mechanism)
    _, point = _foot(mechanism, foot)
    world = RBD.root_frame(mechanism)
    had_contact = false
    contact_before = false
    start_x = 0.0
    max_slip = 0.0
    samples = 0
    for q in qs
        RBD.set_configuration!(state, q)
        p = RBD.transform(state, point, world).v
        contact = p[3] <= ground_tolerance
        if contact
            had_contact = true
            samples += 1
            if contact_before
                max_slip = max(max_slip, abs(p[1] - start_x))
            else
                start_x = p[1]
            end
        end
        contact_before = contact
    end
    slipped = had_contact && max_slip > slip_tolerance
    return (
        was_in_contact = had_contact,
        stayed_in_place = had_contact && !slipped,
        slipped,
        max_slip,
        contact_samples = samples,
    )
end

function _contact_slip_status(
    mechanism,
    qs;
    contact_foot::Symbol = :grounded,
    ground_tolerance = 1e-5,
    slip_tolerance = 1e-4,
)
    contact_foot in (:grounded, :swing, :both) ||
        throw(ArgumentError("contact_foot must be :grounded, :swing, or :both"))
    status(foot) =
        _single_foot_slip_status(mechanism, qs, foot; ground_tolerance, slip_tolerance)
    contact_foot === :grounded && return status(:grounded)
    contact_foot === :swing && return status(:swing)
    return (grounded = status(:grounded), swing = status(:swing))
end

"""Check slip in an already-generated compliant full trajectory."""
function contact_slip_status(model::SimulationModel, qs; kwargs...)
    return _contact_slip_status(_compliant_mechanism(model, 0.0), qs; kwargs...)
end

"""
    simulate_trajectory(model, q0, q̇0, final_time, controller!; foot_model = :locked)

Run either the exact anchored-ankle design model (`:locked`) or the planar
six-DoF compliant-foot model (`:compliant`). The latter adds ankle `x,z` motion
and Hunt--Crossley/Coulomb contact at both feet. `controller!(u, t, q, q̇)`
fills the direct three-motor input `[τ2, τ3, τ4]`.
"""
function simulate_trajectory(
    model::SimulationModel,
    q0,
    q̇0,
    final_time,
    controller!;
    foot_model::Symbol = :locked,
    μ = 0.8,
    Δt = 1e-3,
    normal_stiffness = 1e5,
    normal_damping = 1.0,
    tangential_stiffness = 5e4,
    tangential_damping = 1e4,
    return_contact_status = false,
    contact_foot::Symbol = :grounded,
    ground_tolerance = 1e-5,
    slip_tolerance = 1e-4,
)
    foot_model in (:locked, :compliant) ||
        throw(ArgumentError("foot_model must be :locked or :compliant"))
    mechanism =
        foot_model === :locked ? model.mechanism :
        _compliant_mechanism(
            model,
            μ;
            normal_stiffness,
            normal_damping,
            tangential_stiffness,
            tangential_damping,
        )
    state = state_from_coordinates(model, q0, q̇0; mechanism)
    u = zeros(3)
    function full_controller!(τ, t, current_state)
        q, q̇ = coordinates_from_state(model, current_state)
        controller!(u, t, q, q̇)
        return _set_full_torques!(τ, current_state, model, u)
    end
    ts, qs, vs = RBD.simulate(state, final_time, full_controller!; Δt)
    reduced = [
        begin
            statek = RBD.MechanismState(mechanism)
            RBD.set_configuration!(statek, q)
            RBD.set_velocity!(statek, v)
            coordinates_from_state(model, statek)
        end for (q, v) in zip(qs, vs)
    ]
    qhist = [x[1] for x in reduced]
    q̇hist = [x[2] for x in reduced]
    if return_contact_status
        status =
            foot_model === :compliant ?
            _contact_slip_status(
                mechanism,
                qs;
                contact_foot,
                ground_tolerance,
                slip_tolerance,
            ) : nothing
        return ts, qhist, q̇hist, qs, vs, status
    end
    return ts, qhist, q̇hist, qs, vs
end

"""Create position, velocity, and optional anchored-ankle-force plots."""
function plot_trajectory(ts, qhist, q̇hist; forces = nothing)
    qplot = Plots.plot(
        ts,
        permutedims(reduce(hcat, qhist));
        label = reshape(collect(COORDINATE_NAMES), 1, :),
        xlabel = "time (s)",
        ylabel = "angle (rad)",
    )
    vplot = Plots.plot(
        ts,
        permutedims(reduce(hcat, q̇hist));
        label = reshape(collect(COORDINATE_NAMES), 1, :),
        xlabel = "time (s)",
        ylabel = "angular velocity (rad/s)",
    )
    fplot =
        forces === nothing ? nothing :
        Plots.plot(
            ts,
            permutedims(reduce(hcat, forces));
            label = reshape(["Fx", "Fy", "Fz"], 1, :),
            xlabel = "time (s)",
            ylabel = "force (N)",
        )
    return (positions = qplot, velocities = vplot, forces = fplot)
end

"""Open MeshCat and animate a locked or compliant full trajectory."""
function animate_trajectory(
    model::SimulationModel,
    ts,
    qs;
    foot_model::Symbol = :locked,
    μ = 0.8,
    realtimerate = 1.0,
)
    mechanism = foot_model === :locked ? model.mechanism : _compliant_mechanism(model, μ)
    vis = MeshCatMechanisms.MechanismVisualizer(
        mechanism,
        MeshCatMechanisms.URDFVisuals(model.urdf_file),
    )
    MeshCatMechanisms.animate(vis, ts, qs; realtimerate)
    return vis
end

end # module
