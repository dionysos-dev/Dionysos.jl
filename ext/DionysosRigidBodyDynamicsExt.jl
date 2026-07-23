module DionysosRigidBodyDynamicsExt

import Dionysos
const ST = Dionysos.System

using RigidBodyDynamics
using MeshCat
using MeshCatMechanisms

function Dionysos.animate_mechanism_trajectory(
    urdf::AbstractString,
    x_traj::ST.Trajectory;
    configuration::Function = x -> x,
    joint_names = nothing,
    Δt = 0.05,
    fps = nothing,
    open_visualizer::Bool = true,
    play::Bool = true,
)
    mechanism = parse_urdf(urdf)
    state = MechanismState(mechanism)

    if joint_names === nothing
        joints_to_use = collect(joints(mechanism))
    else
        joints_to_use = map(joint_names) do joint_name
            joint = findjoint(mechanism, joint_name)
            joint === nothing && error("Could not find joint `$joint_name` in URDF.")
            return joint
        end
    end

    xs = collect(ST.states(x_traj))
    isempty(xs) && error("x_traj is empty")

    fps === nothing && (fps = round(Int, 1 / Δt))

    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
    vis = mvis.visualizer

    open_visualizer && open(vis)

    anim = MeshCat.Animation(vis; fps = fps)

    for k in eachindex(xs)
        q = collect(configuration(xs[k]))

        length(q) == length(joints_to_use) || error(
            "configuration(x) returned $(length(q)) values, but $(length(joints_to_use)) joints are used.",
        )

        for (joint, qi) in zip(joints_to_use, q)
            set_configuration!(state, joint, qi)
        end

        MeshCat.atframe(anim, k) do
            return MeshCatMechanisms.set_configuration!(mvis, configuration(state))
        end
    end

    MeshCat.setanimation!(vis, anim; play = play)

    return vis, anim
end

end
