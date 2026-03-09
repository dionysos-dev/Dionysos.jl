import Dionysos

const ST = Dionysos.System

export CandidateTrajectory, horizon, n_states, states, inputs
# pas certain que ça devrait exister
struct CandidateTrajectory{TX<:ST.Trajectory, TU<:ST.Trajectory, M}
    x_traj::TX
    u_traj::TU
    Ts::Float64
    source::Symbol
    metadata::M
end

function CandidateTrajectory(
    x_traj::ST.Trajectory,
    u_traj::ST.Trajectory;
    Ts::Real,
    source::Symbol = :unknown,
    metadata = (;),
)
    Tsf = Float64(Ts)
    Tsf > 0 || error("Ts must be > 0")

    nx = length(x_traj)
    nu = length(u_traj)
    nx == nu + 1 || error("Expected length(x_traj) == length(u_traj) + 1")
    nx >= 2 || error("Trajectory too short")

    return CandidateTrajectory{typeof(x_traj), typeof(u_traj), typeof(metadata)}(
        x_traj,
        u_traj,
        Tsf,
        source,
        metadata,
    )
end

horizon(c::CandidateTrajectory) = length(c.u_traj)
n_states(c::CandidateTrajectory) = length(c.x_traj)
states(c::CandidateTrajectory) = c.x_traj
inputs(c::CandidateTrajectory) = c.u_traj

function Base.show(io::IO, c::CandidateTrajectory)
    print(
        io,
        "CandidateTrajectory(source=",
        c.source,
        ", Ts=",
        c.Ts,
        ", n_states=",
        n_states(c),
        ", horizon=",
        horizon(c),
        ")",
    )
end
