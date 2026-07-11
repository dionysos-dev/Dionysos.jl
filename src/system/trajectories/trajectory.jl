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

Result of a closed-loop simulation: the state trajectory `x`, the input
trajectory `u`, and — for dynamic controllers — the controller-memory
trajectory `q` (`nothing` for static controllers). Destructures as
`x_traj, u_traj[, q_traj] = traj`.
"""
struct ClosedLoopTrajectory{XT, UT, QT}
    x::Trajectory{XT}
    u::Trajectory{UT}
    q::QT
end

ClosedLoopTrajectory(x::Trajectory, u::Trajectory) = ClosedLoopTrajectory(x, u, nothing)

function Base.iterate(traj::ClosedLoopTrajectory, state::Int = 1)
    state == 1 && return (traj.x, 2)
    state == 2 && return (traj.u, 3)
    state == 3 && traj.q !== nothing && return (traj.q, 4)
    return nothing
end

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
