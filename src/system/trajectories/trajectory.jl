export ContinuousTrajectory, ContinuousTrajectoryAttribute
export AutomatonPath
export Trajectory
export last_mode

"""
    AutomatonPath{Q, TT}

A path through a hybrid automaton: a starting mode `q_0` and a sequence of discrete
`transitions`. This is a *search-space object* — the partial candidate grown (via
[`append`](@ref)) by the branch-and-bound solver — not a state/input rollout over time
(that is a `Trajectory`).
"""
struct AutomatonPath{Q, TT}
    q_0::Q
    transitions::Vector{TT}
end

function AutomatonPath{TT}(q_0) where {TT}
    return AutomatonPath(q_0, TT[])
end

function last_mode(system, traj::AutomatonPath)
    if isempty(traj.transitions)
        return traj.q_0
    else
        return HybridSystems.target(system, traj.transitions[end])
    end
end
Base.length(traj::AutomatonPath) = length(traj.transitions)
function append(traj::AutomatonPath, t)
    return AutomatonPath(traj.q_0, [traj.transitions; t])
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
    Trajectory{S, I, T, M, Q}

A system trajectory stored as parallel channels. The `states` channel is always
present; every other channel is optional and is `nothing` when absent:

- `inputs`  — applied inputs; `nothing` for an open-loop / state-only trajectory,
- `times`   — physical time at each step (timed systems),
- `modes`   — active discrete mode at each step (hybrid systems),
- `memory`  — controller memory at each step (dynamic controllers).

So a closed-loop rollout is just a `Trajectory` that carries `inputs`, and a plain
state trajectory is one with `inputs === nothing`. Build it with
`Trajectory(states; inputs = ..., times = ..., modes = ..., memory = ...)` and read
channels with [`states`](@ref), [`inputs`](@ref), [`times`](@ref), [`modes`](@ref),
[`memory`](@ref).
"""
struct Trajectory{S, I, T, M, Q}
    states::Vector{S}
    inputs::I
    times::T
    modes::M
    memory::Q
end

function Trajectory(
    states::AbstractVector;
    inputs = nothing,
    times = nothing,
    modes = nothing,
    memory = nothing,
)
    return Trajectory(states, inputs, times, modes, memory)
end

# Channel accessors — read a single channel as a plain vector (or `nothing` when
# absent). Deliberately unexported: `states`/`modes` would otherwise clash with
# `HybridSystems.states`/`modes` at unqualified call sites.
"State channel (vector of visited states) of a trajectory."
states(traj::Trajectory) = traj.states
"Input channel (vector of applied inputs), or `nothing` for an open-loop trajectory."
inputs(traj::Trajectory) = traj.inputs
"Controller-memory channel, or `nothing` for static controllers."
memory(traj::Trajectory) = traj.memory
"Physical-time channel, or `nothing` when the trajectory is untimed."
times(traj::Trajectory) = traj.times
"Discrete-mode channel, or `nothing` when the trajectory is non-hybrid."
modes(traj::Trajectory) = traj.modes

Base.length(traj::Trajectory) = length(traj.states)

@recipe function f(traj::Trajectory; dims = [1, 2], arrows = true)
    traj_label = get(plotattributes, :label, "")
    xs = states(traj)

    # first point carries the user label, the rest stay unlabelled
    @series begin
        dims := dims
        label := traj_label
        UT.DrawPoint(xs[1])
    end

    for k in 1:(length(xs) - 1)
        @series begin
            dims := dims
            label := ""
            UT.DrawPoint(xs[k + 1])
        end
        if arrows
            @series begin
                dims := dims
                label := ""
                UT.DrawArrow(xs[k], xs[k + 1])
            end
        end
    end
end

"""
    get_cost_trajectory(traj::Trajectory, c) -> (cost, total_cost)

Evaluate the stage cost `c(state, input)` along `traj` (which must carry an
`inputs` channel). Returns the per-step `cost` vector and the summed `total_cost`.
"""
function get_cost_trajectory(traj::Trajectory, c)
    xs = states(traj)
    us = inputs(traj)
    us === nothing && error("get_cost_trajectory needs a trajectory with inputs")
    @assert length(xs) == length(us) + 1

    cs = Float64[]
    total_cost = 0.0

    for i in eachindex(us)
        ci = c(xs[i], us[i])
        total_cost += ci
        push!(cs, ci)
    end

    return (cost = cs, total_cost = total_cost)
end
