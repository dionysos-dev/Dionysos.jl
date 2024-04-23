using BipedRobot
using RigidBodyDynamics

@static if get(ENV, "CI", "false") == "false"
    using MeshCatMechanisms
    using Plots
    using BenchmarkTools
end

urdf = joinpath(dirname(dirname(pathof(BipedRobot))), "deps", "doublependulum.urdf")
doublependulum = parse_urdf(Float64, urdf)

const state = MechanismState(doublependulum)

@static if get(ENV, "CI", "false") == "false"
    vis = MechanismVisualizer(doublependulum, URDFVisuals(urdf))
    #IJuliaCell(vis)
end

set_configuration!(state, [1.0, -1.5])
@static if get(ENV, "CI", "false") == "false"
    set_configuration!(vis, configuration(state))
end

ts, qs, vs = simulate(state, 5.0; Δt = 1e-3);

@static if get(ENV, "CI", "false") == "false"
    MeshCatMechanisms.animate(vis, ts, qs)
end

shoulder_angles = collect(q[1] for q in qs)

@static if get(ENV, "CI", "false") == "false"
    plot(
        ts,
        shoulder_angles;
        xlabel = "Time [s]",
        ylabel = "Angle [rad]",
        label = "Shoulder",
    )
end

bodies(doublependulum)
joints(doublependulum)
fixedjoint, shoulder, elbow = joints(doublependulum)
@show fixedjoint shoulder elbow;
q = configuration(state)
v = velocity(state)
@show q v;
q[1]
q[shoulder]

frame_after(elbow)
@static if get(ENV, "CI", "false") == "false"
    setelement!(vis, frame_after(elbow))
end
frame_before(shoulder)
@static if get(ENV, "CI", "false") == "false"
    setelement!(vis, frame_before(shoulder))
end
p = Point3D(frame_after(elbow), 0.0, 0.0, -2.0)

radius = 0.1
@static if get(ENV, "CI", "false") == "false"
    setelement!(vis, p, radius, "tip")
end
p = transform(state, p, root_frame(doublependulum))

displacement = FreeVector3D(frame_after(elbow), 2.0, 3.0, 4.0)

try
    p + displacement
catch err
    println("failed!")
    err
end

T = transform_to_root(state, frame_after(elbow))
relative_transform(state, frame_after(elbow), frame_after(shoulder))

# center_of_mass(state)

twist = relative_twist(state, frame_after(elbow), frame_before(shoulder))

transform(state, twist, frame_after(elbow))

M = mass_matrix(state)

mass_matrix!(M, state)

@static if get(ENV, "CI", "false") == "false"
    @btime (setdirty!($state); mass_matrix!($M, $state))
else
    mass_matrix!(M, state)
end

v̇ = similar(velocity(state))
v̇ .= [2.0; 3.0] # the joint acceleration vector, i.e., the time derivative of the joint velocity vector v
τ = inverse_dynamics(state, v̇)
@show τ;

result = DynamicsResult{Float64}(doublependulum);

dynamics!(result, state)
@show result.v̇;

# *********************************************************************************************************************

using Symbolics
@variables q1 q2 v1 v2 real = true
q = [q1, q2]
v = [v1, v2]

simplify.(mass_matrix(MechanismState(doublependulum, q, v)))
