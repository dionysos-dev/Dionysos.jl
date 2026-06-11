import Dionysos

const ST = Dionysos.System

export CandidateTrajectory, horizon, n_states, states, inputs

struct CandidateTrajectory{TX <: ST.Trajectory, TU <: ST.Trajectory, M}
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
    @assert Tsf > 0.0

    nx = length(x_traj)
    nu = length(u_traj)
    @assert nx == nu + 1
    @assert nx >= 2

    return CandidateTrajectory(x_traj, u_traj, Tsf, source, metadata)
end

horizon(c::CandidateTrajectory) = length(c.u_traj)
n_states(c::CandidateTrajectory) = length(c.x_traj)
states(c::CandidateTrajectory) = c.x_traj
inputs(c::CandidateTrajectory) = c.u_traj

function Base.show(io::IO, c::CandidateTrajectory)
    return print(
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
