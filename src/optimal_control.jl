export ZeroCost, ConstantCost, QuadraticControlCost
export OptimalControlProblem

struct DiscreteTrajectory{TT}
    q_0::Int
    transitions::Vector{TT}
end
function DiscreteTrajectory{TT}(q_0::Int) where TT
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
    DiscreteTrajectory(traj.q_0, [traj.transitions; t])
end

struct ContinuousTrajectory{T, VT<:AbstractVector{T}}
    x::Vector{VT}
    u::Vector{VT}
end

struct HybridTrajectory{T, TT, VT<:AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T, VT}
end

struct ZeroCost end

struct ConstantCost{T}
    cost::T
end

struct QuadraticControlCost{T, MT <: AbstractMatrix{T}}
    Q::MT
end

struct OptimalControlProblem{S, Q, X0, XC, TC}
    system::S
    q_0::Q
    x_0::X0
    state_cost::XC
    transition_cost::TC
    q_T::Q
    number_of_time_steps::Int
end

export optimal_control
function optimal_control end

function _zero_steps(prob)
    if prob.q_T == prob.q_0
        return ContinuousTrajectory(Vector{Float64}[], Vector{Float64}[]), 0.0
    else
        return
    end
end
