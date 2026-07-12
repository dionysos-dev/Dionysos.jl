export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export last_mode

"""
    DiscreteTrajectory{Q, TT}

`q_0` is the starting mode and `transitions` is a sequence of discrete transitions in the system.
"""
struct DiscreteTrajectory{Q, TT}
    q_0::Q
    transitions::Vector{TT}
end

function DiscreteTrajectory{TT}(q_0) where {TT}
    return DiscreteTrajectory(q_0, TT[])
end

function last_mode(system, traj::DiscreteTrajectory)
    if isempty(traj.transitions)
        return traj.q_0
    else
        return HybridSystems.target(system, traj.transitions[end])
    end
end
Base.length(traj::DiscreteTrajectory) = length(traj.transitions)
function append(traj::DiscreteTrajectory, t)
    return DiscreteTrajectory(traj.q_0, [traj.transitions; t])
end

"""
    ContinuousTrajectory{T, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}

`x` is a sequence of points in the state space and `u` is a sequence of points in the input space.
"""
struct ContinuousTrajectory{T, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    x::Vector{XVT}
    u::Vector{UVT}
end

struct ContinuousTrajectoryAttribute <: MOI.AbstractModelAttribute end

"""
    Trajectory{T}

provides the sequence of some elements of a trajectory.
"""
struct Trajectory{T}
    seq::Vector{T}
end

"""
    ClosedLoopTrajectory

Result of a closed-loop simulation, stored as parallel channels. `x` is the
state trajectory and `u` the input trajectory (always present); `q` is the
controller-memory trajectory for dynamic controllers (`nothing` for static
ones). The optional `times` and `modes` channels carry, respectively, the
physical time and the discrete mode at each step for timed / hybrid rollouts
(`nothing` when absent). Destructures as `x_traj, u_traj[, q_traj] = traj`;
read individual channels with `states`, `inputs`, `memory`, `times`, `modes`.
"""
struct ClosedLoopTrajectory{XT, UT, QT, TT, MT}
    x::Trajectory{XT}
    u::Trajectory{UT}
    q::QT
    times::TT
    modes::MT
end

function ClosedLoopTrajectory(
    x::Trajectory,
    u::Trajectory,
    q = nothing;
    times = nothing,
    modes = nothing,
)
    return ClosedLoopTrajectory(x, u, q, times, modes)
end

function Base.iterate(traj::ClosedLoopTrajectory, state::Int = 1)
    state == 1 && return (traj.x, 2)
    state == 2 && return (traj.u, 3)
    state == 3 && traj.q !== nothing && return (traj.q, 4)
    return nothing
end

# Channel accessors — read a single channel as a plain vector (or `nothing`
# when absent). Deliberately unexported: `states`/`modes` would otherwise clash
# with `HybridSystems.states`/`modes` at unqualified call sites.
"State channel (vector of visited states) of a closed-loop trajectory."
states(traj::ClosedLoopTrajectory) = traj.x.seq
"Input channel (vector of applied inputs) of a closed-loop trajectory."
inputs(traj::ClosedLoopTrajectory) = traj.u.seq
"Controller-memory channel, or `nothing` for static controllers."
memory(traj::ClosedLoopTrajectory) = traj.q === nothing ? nothing : traj.q.seq
"Physical-time channel, or `nothing` when the trajectory is untimed."
times(traj::ClosedLoopTrajectory) = traj.times
"Discrete-mode channel, or `nothing` when the trajectory is non-hybrid."
modes(traj::ClosedLoopTrajectory) = traj.modes

Base.length(traj::Trajectory) = length(traj.seq)
get_elem(traj::Trajectory, n::Int) = traj.seq[n]
enum_elems(traj::Trajectory) = traj.seq

@recipe function f(traj::Trajectory; dims = [1, 2], arrows = true)
    traj_label = get(plotattributes, :label, "")

    # first point carries the user label, the rest stay unlabelled
    @series begin
        dims := dims
        label := traj_label
        UT.DrawPoint(traj.seq[1])
    end

    for k in 1:(length(traj.seq) - 1)
        @series begin
            dims := dims
            label := ""
            UT.DrawPoint(traj.seq[k + 1])
        end
        if arrows
            @series begin
                dims := dims
                label := ""
                UT.DrawArrow(traj.seq[k], traj.seq[k + 1])
            end
        end
    end
end

function get_cost_trajectory(x_traj::Trajectory, u_traj::Trajectory, c)
    @assert length(x_traj) == length(u_traj) + 1

    cs = Float64[]
    total_cost = 0.0

    for i in 1:length(u_traj)
        ci = c(x_traj.seq[i], u_traj.seq[i])
        total_cost += ci
        push!(cs, ci)
    end

    return (c = Trajectory(cs), total_cost)
end
