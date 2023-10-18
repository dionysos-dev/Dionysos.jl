export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export last_mode

using JuMP
using HybridSystems
using Polyhedra

"""
    DiscreteTrajectory{Q, TT}

`q_0` is the starting mode and `transitions` is a sequence of discrete
transitions in the system
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
        return target(system, traj.transitions[end])
    end
end
Base.length(traj::DiscreteTrajectory) = length(traj.transitions)
function append(traj::DiscreteTrajectory, t)
    return DiscreteTrajectory(traj.q_0, [traj.transitions; t])
end

"""
    ContinuousTrajectory{T, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}

`x` is a sequence of points in the state space and `u` is a sequence of points
in the input space
"""
struct ContinuousTrajectory{T, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    x::Vector{XVT}
    u::Vector{UVT}
end

struct ContinuousTrajectoryAttribute <: MOI.AbstractModelAttribute end

struct HybridTrajectory{T, TT, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T, XVT, UVT}
end

# Trajectory closed loop
function get_closed_loop_trajectory(contsys, controller, x0, nstep; stopping = (x) -> false)
    x = x0
    x_traj = [x]
    u_traj = []
    i = 0
    while !stopping(x) && i ≤ nstep
        u = controller(x)
        if u === nothing
            break
        end
        x = contsys.sys_map(x, u, contsys.tstep)
        push!(x_traj, x)
        push!(u_traj, u)
        i = i + 1
    end
    return x_traj, u_traj
end

function get_closed_loop_trajectory(
    f_eval,
    c_eval,
    cost_eval,
    x0,
    nstep;
    stopping = (x) -> false,
    noise = false,
)
    x = x0
    x_traj = [x]
    u_traj = []
    cost_traj = []
    i = 0
    while !stopping(x) && i ≤ nstep
        u = c_eval(x)
        if u === nothing
            break
        end
        push!(u_traj, u)
        push!(cost_traj, cost_eval(x, u))
        if noise
            w = zeros(2)
            x = f_eval(x, u, w)
        else
            x = f_eval(x, u)
        end
        push!(x_traj, x)
        i = i + 1
    end
    return x_traj, u_traj, cost_traj
end
