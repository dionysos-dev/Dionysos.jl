export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export last_mode

using JuMP
using HybridSystems
using Polyhedra

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
        return target(system, traj.transitions[end])
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
    HybridTrajectory{T, TT, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}

`discrete` is the discrete trajectory of type `DiscreteTrajectory` and  `continuous` is a the continuous trajectory of type `ContinuousTrajectory`.
"""
struct HybridTrajectory{T, TT, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T, XVT, UVT}
end

"""
    Trajectory{T}

provides the sequence of some elements of a trajectory.
"""
struct Trajectory{T}
    seq::Vector{T}
end

Base.length(traj::Trajectory) = length(traj.seq)
get_elem(traj::Trajectory, n::Int) = traj.seq[n]

@recipe function f(traj::Trajectory)
    return UT.DrawTrajectory(traj.seq)
end

"""
    Control_trajectory{T1, T2}

provides the sequence of states and inputs of a trajectory.
"""
struct Control_trajectory{T1, T2}
    states::Trajectory{T1}
    inputs::Trajectory{T2}
end

Base.length(traj::Control_trajectory) = length(traj.states)
get_state(traj::Control_trajectory, n::Int) = get_elem(traj.states, n)
get_input(traj::Control_trajectory, n::Int) = get_elem(traj.inputs, n)
get_elem(traj::Control_trajectory, n::Int) = (get_state(traj, n), get_input(traj, n))

@recipe function f(traj::Control_trajectory)
    return traj.states
end

"""
    Cost_control_trajectory{T1, T2, T3}

provides the sequence of states, inputs (via `Control_trajectory`) and costs of a trajectory.
"""
struct Cost_control_trajectory{T1, T2, T3}
    control_trajectory::Control_trajectory{T1, T2}
    costs::Trajectory{T3}
end

Base.length(traj::Cost_control_trajectory) = length(traj.control_trajectory)
get_state(traj::Cost_control_trajectory, n::Int) = get_state(traj.control_trajectory, n)
get_input(traj::Cost_control_trajectory, n::Int) = get_input(traj.control_trajectory, n)
get_cost(traj::Cost_control_trajectory, n::Int) = get_elem(traj.costs, n)
get_elem(traj::Cost_control_trajectory, n::Int) =
    (get_state(traj, n), get_input(traj, n), get_cost(traj, n))

@recipe function f(traj::Cost_control_trajectory)
    return traj.control_trajectory
end

function get_closed_loop_trajectory(contsys, controller, x0, nstep; stopping = (x) -> false)
    x = x0
    x_traj = [x]
    u_traj = []
    u_traj2 = []
    t_traj = []
    i = 0
    total_time = 0.0
    energy     = 0.0
    control_effort = 0.0    
    while !stopping(x) && i ≤ nstep
        u, p = controller(x) # u, tstep_cur = controller(x)
        timestep = contsys.tstep + p * 0.05 #contsys.tstep * 1.5 ^ p
        println("u = $u, p = $p, timestep = $timestep")
        energy += x[2] * u[1] * timestep
        control_effort += u[1]^2 * timestep
        if u === nothing
            break
        end
        push!(t_traj, total_time)
        total_time += timestep
        x = contsys.sys_map(x, u, timestep)
        push!(x_traj, x)
        push!(u_traj, u)
        push!(u_traj2, u[1])
        i = i + 1
    end
    # path = "C:/Users/adrie/OneDrive - UCL/Master 2/mémoire visus/data/"
    # open(joinpath(path, "u.txt"), "w") do file
    #     for i in 1:length(u_traj2)
    #         t = t_traj[i]
    #         u = u_traj2[i]
    #         println(file, "$t $u")
    #     end
    # end
    println("Total time: $total_time")
    println("Energy: $energy")
    println("Control effort: $control_effort")
    return Control_trajectory(Trajectory(x_traj), Trajectory(u_traj))
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
    control_trajectory = Control_trajectory(Trajectory(x_traj), Trajectory(u_traj))
    return Cost_control_trajectory(control_trajectory, Trajectory(cost_traj))
end
