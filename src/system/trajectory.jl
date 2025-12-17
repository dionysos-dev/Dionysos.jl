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

@recipe function f(traj::Trajectory; dims=[1,2])
    @series begin
        dims := dims
        UT.DrawTrajectory(traj.seq)
    end
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

@recipe function f(traj::Control_trajectory; dims=[1,2])
    @series begin
        dims := dims
        traj.states
    end
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

function get_cost(traj::Cost_control_trajectory)
    isempty(traj.costs.seq) && return 0.0
    return sum(traj.costs.seq)
end

@recipe function f(traj::Cost_control_trajectory; dims=[1,2])
    @series begin
        dims := dims
        traj.control_trajectory
    end
end

"""
    wrap_coord(x::SVector{N, T}, periodic_dims::SVector{P, Int}, periods::SVector{P, T}; start = zeros(SVector{P, T}))

Wraps the vector `x` into a periodic domain along dimensions specified in `periodic_dims`,
with period lengths `periods` and optional offset `start`.

# Arguments
- `x`: The coordinate vector to wrap.
- `periodic_dims`: Indices of the periodic dimensions.
- `periods`: Period lengths for the periodic dimensions.
- `start` (optional): Starting values of the periodic domains (defaults to `0.0`).

# Returns
A wrapped `SVector` where each periodic dimension is mapped to the interval `[start[i], start[i] + periods[i])`.
"""
function wrap_coord(
    x::SVector{N, T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {N, P, T}
    return SVector{N, T}(ntuple(d -> begin
        i = findfirst(isequal(d), periodic_dims)
        if i === nothing
            x[d]
        else
            s = start[i]
            p = periods[i]
            mod(x[d] - s, p) + s
        end
    end, N))
end

function get_periodic_wrapper(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {P, T}
    return x -> wrap_coord(x, periodic_dims, periods; start = start)
end

function get_closed_loop_trajectory(
    system::MS.ConstrainedBlackBoxControlDiscreteSystem,
    controller::ContinuousController,
    x0,
    nstep;
    stopping = (x) -> false,
    periodic_wrapper = x -> x,
)
    x = periodic_wrapper(x0)
    x_traj, u_traj = [x], []

    for _ in 1:nstep
        stopping(x) && break
        u = get_control(controller, x)
        u === nothing && break
        x = MS.mapping(system)(x, u)
        x = periodic_wrapper(x)
        push!(x_traj, x)
        push!(u_traj, u)
    end

    return Control_trajectory(Trajectory(x_traj), Trajectory(u_traj))
end

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
    return Control_trajectory(Trajectory(x_traj), Trajectory(u_traj))
end

function get_closed_loop_trajectory(
    f_eval,
    controller::ContinuousController,
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
        u = get_control(controller, x)
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
