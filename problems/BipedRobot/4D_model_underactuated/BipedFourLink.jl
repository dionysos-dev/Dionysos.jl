"""Planar four-coordinate biped: exact ankle support for design and soft feet for simulation."""
module BipedFourLink

using LinearAlgebra
using StaticArrays
using RigidBodyDynamics
using Plots
using MeshCatMechanisms

const RBD = RigidBodyDynamics

export BipedModel, build_model, state_from_coordinates!, state_from_coordinates,
    coordinates_from_state, world_position, foot_position,
    coordinate_mass_matrix, bias_forces, accelerations, coordinate_inverse_dynamics,
    velocity_controlled_derivative, velocity_controlled_contact_estimate,
    ankle_reaction_wrench, ankle_reaction_force, in_friction_cone,
    simulate_trajectory, contact_slip_status, plot_trajectory, animate_trajectory,
    COORDINATE_NAMES, INPUT_NAMES

"""
Four-coordinate design model.

`mechanism` has an exact revolute ankle support; `q1_joint` is its tree joint.
The compliant simulation model is constructed separately with a planar root.
"""
struct BipedModel
    mechanism
    q1_joint
    joint_names::NTuple{3, String}
    urdf_file::String
end

const COORDINATE_NAMES = ("q1_ankle", "q2_grounded_knee", "q3_hip", "q4_swing_knee")
const INPUT_NAMES = ("τ2_grounded_knee", "τ3_hip", "τ4_swing_knee")

_parse(urdf_file, root_joint_type) = RBD.URDF.parse_urdf(
    urdf_file;
    root_joint_type,
    remove_fixed_tree_joints = false,
    gravity = SVector(0.0, 0.0, -9.81),
)

"""
    build_model(urdf_file)

Create the four-DoF design mechanism. `q1` is a ground-ankle revolute joint;
`q2`, `q3`, and `q4` are the grounded knee, hip, and swing knee.
"""
function build_model(urdf_file::AbstractString)
    mechanism = _parse(urdf_file, RBD.Revolute(SVector(0.0, 1.0, 0.0)))
    BipedModel(
        mechanism,
        RBD.findjoint(mechanism, "grounded_tibia_to_world"),
        ("grounded_knee", "hip", "swing_knee"),
        String(urdf_file),
    )
end

_joint(model, index, mechanism = model.mechanism) = RBD.findjoint(mechanism, model.joint_names[index])

"""Set a locked-ankle state from `q = [q1, q2, q3, q4]` and `q̇`."""
function state_from_coordinates!(state, model::BipedModel, q::AbstractVector, q̇::AbstractVector = zeros(4))
    length(q) == 4 || throw(DimensionMismatch("q must have four coordinates"))
    length(q̇) == 4 || throw(DimensionMismatch("q̇ must have four coordinates"))
    RBD.set_configuration!(state, model.q1_joint, q[1])
    RBD.set_velocity!(state, model.q1_joint, q̇[1])
    for i in 1:3
        RBD.set_configuration!(state, _joint(model, i), q[i + 1])
        RBD.set_velocity!(state, _joint(model, i), q̇[i + 1])
    end
    state
end

function state_from_coordinates(model::BipedModel, q::AbstractVector, q̇::AbstractVector = zeros(4))
    state = RBD.MechanismState(model.mechanism)
    state_from_coordinates!(state, model, q, q̇)
end

"""Extract `[q1, q2, q3, q4]` and its derivative from either internal model."""
function coordinates_from_state(model::BipedModel, state)
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
    qs = [q1; [q[RBD.configuration_range(state, _joint(model, i, mechanism))][1] for i in 1:3]]
    q̇s = [q̇1; [v[RBD.velocity_range(state, _joint(model, i, mechanism))[1]] for i in 1:3]]
    qs, q̇s
end

"""World location of a URDF joint or body."""
function world_position(model::BipedModel, state, name::AbstractString)
    mechanism = state.mechanism
    frame = try
        RBD.frame_after(RBD.findjoint(mechanism, name))
    catch
        RBD.default_frame(RBD.findbody(mechanism, name))
    end
    RBD.translation(RBD.transform_to_root(state, frame))
end

"""World position of the grounded or swing foot point (`foot = :grounded` or `:swing`)."""
function foot_position(model::BipedModel, state; foot::Symbol = :grounded)
    foot in (:grounded, :swing) || throw(ArgumentError("foot must be :grounded or :swing"))
    mechanism = state.mechanism
    body_name, local_z = foot === :grounded ? ("grounded_tibia", 0.0) : ("swing_tibia", -0.172)
    body = RBD.findbody(mechanism, body_name)
    point = RBD.Point3D(RBD.default_frame(body), 0.0, 0.0, local_z)
    RBD.transform(state, point, RBD.root_frame(mechanism)).v
end

function _full_torques(state, model, u)
    length(u) == 3 || throw(DimensionMismatch("u must be [τ2, τ3, τ4]"))
    τ = zeros(eltype(RBD.configuration(state)), RBD.num_velocities(state.mechanism))
    _set_full_torques!(τ, state, model, u)
end

function _set_full_torques!(τ, state, model, u)
    length(u) == 3 || throw(DimensionMismatch("u must be [τ2, τ3, τ4]"))
    fill!(τ, zero(eltype(τ)))
    for i in 1:3
        τ[RBD.velocity_range(state, _joint(model, i, state.mechanism))[1]] = u[i]
    end
    τ
end

function _terms(state, model)
    result = RBD.DynamicsResult(state.mechanism)
    RBD.dynamics!(result, state, zeros(RBD.num_velocities(state.mechanism)))
    Matrix(result.massmatrix), collect(parent(result.dynamicsbias)), result
end

"""Return the 4 × 4 mass matrix in `[q1, q2, q3, q4]` coordinates."""
function coordinate_mass_matrix(model::BipedModel, q, q̇ = zeros(4))
    state = state_from_coordinates(model, q, q̇)
    M, _, _ = _terms(state, model)
    Symmetric(M)
end

"""Return gravity and velocity-dependent generalized forces in four coordinates."""
function bias_forces(model::BipedModel, q, q̇ = zeros(4))
    state = state_from_coordinates(model, q, q̇)
    _, c, _ = _terms(state, model)
    c
end

"""Compute `q̈` under direct motor torques `u = [τ2, τ3, τ4]`; q1 is unactuated."""
function accelerations(model::BipedModel, q, q̇, u)
    state = state_from_coordinates(model, q, q̇)
    result = RBD.DynamicsResult(model.mechanism)
    RBD.dynamics!(result, state, _full_torques(state, model, u))
    collect(result.v̇)
end

"""Return `(u, ankle_torque)` required to realize the prescribed `q̈`."""
function coordinate_inverse_dynamics(model::BipedModel, q, q̇, q̈)
    length(q̈) == 4 || throw(DimensionMismatch("q̈ must have four elements"))
    M, c, _ = _terms(state_from_coordinates(model, q, q̇), model)
    required = M * q̈ + c
    (u = required[2:4], ankle_torque = required[1])
end

"""
    velocity_controlled_derivative(model, x, v; input_accelerations = zeros(3))

Return the derivative of the velocity-controlled design state
`x = [q̇₁, q₁, q₂, q₃, q₄]` for commanded actuated-joint velocities
`v = [v₂, v₃, v₄]`.

The resulting equations are `q̇₂:q̇₄ = v₂:v₄`, `q̇₁ = x[1]`, and the exact
unactuated first Euler--Lagrange row
`M₁₁*q̈₁ + M₁a*v̇ + c₁ = 0`. By default `v̇ = 0`, which is the piecewise-
constant velocity command used by the finite-graph design model.
"""
function velocity_controlled_derivative(model::BipedModel, x, v;
        input_accelerations = zeros(3))
    length(x) == 5 || throw(DimensionMismatch("x must be [q̇1, q1, q2, q3, q4]"))
    length(v) == 3 || throw(DimensionMismatch("v must be [v2, v3, v4]"))
    length(input_accelerations) == 3 || throw(DimensionMismatch("input_accelerations must have three elements"))
    q = collect(x[2:5])
    q̇ = [x[1]; v]
    M = coordinate_mass_matrix(model, q, q̇)
    c = bias_forces(model, q, q̇)
    q̈1 = -(dot(view(M, 1, 2:4), input_accelerations) + c[1]) / M[1, 1]
    [q̈1, x[1], v[1], v[2], v[3]]
end

"""
    velocity_controlled_contact_estimate(model, x, v; input_accelerations = zeros(3))

Return an a-posteriori anchored-ankle wrench estimate for the velocity-
controlled design state. The returned named tuple includes the inferred full
acceleration and the motor torques that would realize it in the full dynamics.
"""
function velocity_controlled_contact_estimate(model::BipedModel, x, v;
        input_accelerations = zeros(3))
    ẋ = velocity_controlled_derivative(model, x, v; input_accelerations)
    q = collect(x[2:5])
    q̇ = [x[1]; v]
    q̈ = [ẋ[1]; input_accelerations]
    (
        wrench = ankle_reaction_wrench(model, q, q̇, q̈),
        force = ankle_reaction_force(model, q, q̇, q̈),
        q̈,
        required = coordinate_inverse_dynamics(model, q, q̇, q̈),
    )
end

"""World-frame wrench applied by the world at the grounded ankle for `(q, q̇, q̈)`."""
function ankle_reaction_wrench(model::BipedModel, q, q̇, q̈)
    state = state_from_coordinates(model, q, q̇)
    result = RBD.DynamicsResult(model.mechanism)
    copyto!(result.v̇, q̈)
    RBD.inverse_dynamics!(result.dynamicsbias, result.jointwrenches, result.accelerations, state, result.v̇)
    grounded_tibia = RBD.findbody(model.mechanism, "grounded_tibia")
    wrench = RBD.joint_wrench(result, grounded_tibia)
    RBD.transform(state, wrench, RBD.root_frame(model.mechanism))
end

"""World-frame `[Fx, Fy, Fz]` grounded-ankle reaction force."""
ankle_reaction_force(model::BipedModel, q, q̇, q̈) = collect(RBD.linear(ankle_reaction_wrench(model, q, q̇, q̈)))

"""Test `||F_t|| ≤ μ F_n` and `F_n ≥ 0` for a world-frame contact force."""
function in_friction_cone(force, μ; normal_axis = 3, normal_sign = 1)
    fn = normal_sign * force[normal_axis]
    ft2 = sum(force[i]^2 for i in eachindex(force) if i != normal_axis)
    fn >= 0 && ft2 <= (μ * fn)^2
end

function _compliant_mechanism(model::BipedModel, μ;
        normal_stiffness = 1e5, normal_damping = 1.0,
        tangential_stiffness = 5e4, tangential_damping = 1e4)
    # +x and -z translations give a +y pitch axis: x × (-z) = +y.
    mechanism = _parse(model.urdf_file,
        RBD.Planar(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, -1.0)))
    world = RBD.root_frame(mechanism)
    RBD.add_environment_primitive!(mechanism, RBD.HalfSpace3D(
        RBD.Point3D(world, 0.0, 0.0, 0.0), RBD.FreeVector3D(world, 0.0, 0.0, 1.0)))
    normal = RBD.hunt_crossley_hertz(k = normal_stiffness; α = normal_damping)
    friction = RBD.ViscoelasticCoulombModel(μ, tangential_stiffness, tangential_damping)
    contact = RBD.SoftContactModel(normal, friction)
    for (name, z) in (("grounded_tibia", 0.0), ("swing_tibia", -0.172))
        body = RBD.findbody(mechanism, name)
        RBD.add_contact_point!(body, RBD.ContactPoint(RBD.Point3D(RBD.default_frame(body), 0.0, 0.0, z), contact))
    end
    mechanism
end

function _compliant_initial_state(model, q, q̇, mechanism)
    state = RBD.MechanismState(mechanism)
    root = RBD.findjoint(mechanism, "grounded_tibia_to_world")
    # Base coordinates [x, -z, q1]; grounded ankle starts at the origin.
    RBD.set_configuration!(state, root, [0.0, 0.0, q[1]])
    RBD.set_velocity!(state, root, [0.0, 0.0, q̇[1]])
    for i in 1:3
        RBD.set_configuration!(state, _joint(model, i, mechanism), q[i + 1])
        RBD.set_velocity!(state, _joint(model, i, mechanism), q̇[i + 1])
    end
    state
end

function _single_foot_slip_status(mechanism, qs, foot;
        ground_tolerance = 1e-5, slip_tolerance = 1e-4)
    state = RBD.MechanismState(mechanism)
    body_name, local_z = foot === :grounded ? ("grounded_tibia", 0.0) : ("swing_tibia", -0.172)
    body = RBD.findbody(mechanism, body_name)
    point = RBD.Point3D(RBD.default_frame(body), 0.0, 0.0, local_z)
    world = RBD.root_frame(mechanism)
    had_contact = false; contact_before = false; start_x = 0.0; max_slip = 0.0; samples = 0
    for q in qs
        RBD.set_configuration!(state, q)
        p = RBD.transform(state, point, world).v
        contact = p[3] <= ground_tolerance
        if contact
            had_contact = true; samples += 1
            if contact_before
                max_slip = max(max_slip, abs(p[1] - start_x))
            else
                start_x = p[1]
            end
        end
        contact_before = contact
    end
    slipped = had_contact && max_slip > slip_tolerance
    (was_in_contact = had_contact, stayed_in_place = had_contact && !slipped,
        slipped, max_slip, contact_samples = samples)
end

function _contact_slip_status(mechanism, qs; contact_foot::Symbol = :grounded,
        ground_tolerance = 1e-5, slip_tolerance = 1e-4)
    contact_foot in (:grounded, :swing, :both) || throw(ArgumentError("contact_foot must be :grounded, :swing, or :both"))
    status(foot) = _single_foot_slip_status(mechanism, qs, foot; ground_tolerance, slip_tolerance)
    contact_foot === :grounded && return status(:grounded)
    contact_foot === :swing && return status(:swing)
    (grounded = status(:grounded), swing = status(:swing))
end

"""Check slip in an already-generated compliant full trajectory."""
function contact_slip_status(model::BipedModel, qs; kwargs...)
    _contact_slip_status(_compliant_mechanism(model, 0.0), qs; kwargs...)
end

"""
    simulate_trajectory(model, q0, q̇0, final_time, controller!; foot_model = :locked)

Run either the exact anchored-ankle design model (`:locked`) or the planar
six-DoF compliant-foot model (`:compliant`). The latter adds ankle `x,z` motion
and Hunt--Crossley/Coulomb contact at both feet. `controller!(u, t, q, q̇)`
fills the direct three-motor input `[τ2, τ3, τ4]`.
"""
function simulate_trajectory(model::BipedModel, q0, q̇0, final_time, controller!;
        foot_model::Symbol = :locked, μ = 0.8, Δt = 1e-3,
        normal_stiffness = 1e5, normal_damping = 1.0,
        tangential_stiffness = 5e4, tangential_damping = 1e4,
        return_contact_status = false, contact_foot::Symbol = :grounded,
        ground_tolerance = 1e-5, slip_tolerance = 1e-4)
    foot_model in (:locked, :compliant) || throw(ArgumentError("foot_model must be :locked or :compliant"))
    mechanism = foot_model === :locked ? model.mechanism : _compliant_mechanism(model, μ;
        normal_stiffness, normal_damping, tangential_stiffness, tangential_damping)
    state = foot_model === :locked ? state_from_coordinates(model, q0, q̇0) :
        _compliant_initial_state(model, q0, q̇0, mechanism)
    u = zeros(3)
    function full_controller!(τ, t, current_state)
        q, q̇ = coordinates_from_state(model, current_state)
        controller!(u, t, q, q̇)
        _set_full_torques!(τ, current_state, model, u)
    end
    ts, qs, vs = RBD.simulate(state, final_time, full_controller!; Δt)
    reduced = [begin
        statek = RBD.MechanismState(mechanism)
        RBD.set_configuration!(statek, q); RBD.set_velocity!(statek, v)
        coordinates_from_state(model, statek)
    end for (q, v) in zip(qs, vs)]
    qhist = [x[1] for x in reduced]
    q̇hist = [x[2] for x in reduced]
    if return_contact_status
        status = foot_model === :compliant ? _contact_slip_status(mechanism, qs;
            contact_foot, ground_tolerance, slip_tolerance) : nothing
        return ts, qhist, q̇hist, qs, vs, status
    end
    ts, qhist, q̇hist, qs, vs
end

"""Create position, velocity, and optional anchored-ankle-force plots."""
function plot_trajectory(ts, qhist, q̇hist; forces = nothing)
    qplot = Plots.plot(ts, permutedims(reduce(hcat, qhist));
        label = reshape(collect(COORDINATE_NAMES), 1, :), xlabel = "time (s)", ylabel = "angle (rad)")
    vplot = Plots.plot(ts, permutedims(reduce(hcat, q̇hist));
        label = reshape(collect(COORDINATE_NAMES), 1, :), xlabel = "time (s)", ylabel = "angular velocity (rad/s)")
    fplot = forces === nothing ? nothing : Plots.plot(ts, permutedims(reduce(hcat, forces));
        label = reshape(["Fx", "Fy", "Fz"], 1, :), xlabel = "time (s)", ylabel = "force (N)")
    (positions = qplot, velocities = vplot, forces = fplot)
end

"""Open MeshCat and animate a locked or compliant full trajectory."""
function animate_trajectory(model::BipedModel, ts, qs;
        foot_model::Symbol = :locked, μ = 0.8, realtimerate = 1.0)
    mechanism = foot_model === :locked ? model.mechanism : _compliant_mechanism(model, μ)
    vis = MeshCatMechanisms.MechanismVisualizer(mechanism, MeshCatMechanisms.URDFVisuals(model.urdf_file))
    MeshCatMechanisms.animate(vis, ts, qs; realtimerate)
    vis
end

end # module
