using RigidBodyDynamics

srcdir = dirname(pathof(RigidBodyDynamics))
urdf = joinpath("docs", "src", "examples", "solvers", "Acrobot.urdf")
mechanism = parse_urdf(urdf)

shoulder, elbow = joints(mechanism)
function simple_control!(torques::AbstractVector, t, state::MechanismState)
    torques[velocity_range(state, shoulder)] .= -1 .* velocity(state, shoulder)
    torques[velocity_range(state, elbow)] .= 10 * sin(t)
end;

state = MechanismState(mechanism)
zero_velocity!(state)
set_configuration!(state, shoulder, 0.7)
set_configuration!(state, elbow, -0.8);


final_time = 4.
ts, qs, vs = simulate(state, final_time, simple_control!; Δt = 0.01);
# println(ts)
println(qs)
println(typeof(qs))

# println(qs)

using MeshCatMechanisms

mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf));
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.);
