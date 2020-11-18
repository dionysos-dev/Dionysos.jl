export ZeroFunction, ConstantFunction, QuadraticControlFunction, PolyhedralFunction
export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export OptimalControlProblem
export last_mode, function_value

using JuMP
using HybridSystems

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

struct ContinuousTrajectoryAttribute <: MOI.AbstractModelAttribute
end

struct HybridTrajectory{T, TT, VT<:AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T, VT}
end

struct ZeroFunction end

struct ConstantFunction{T}
    value::T
end
function_value(f::ConstantFunction, x) = f.value
function Base.:+(f::ConstantFunction, g::ConstantFunction)
    return ConstantFunction(f.value + g.value)
end

struct QuadraticControlFunction{T, MT <: AbstractMatrix{T}}
    Q::MT
end

struct AffineFunction{T}
    a::Vector{T}
    β::T
end
function function_value(f::AffineFunction, x)
    return f.a ⋅ x + f.β
end

struct PolyhedralFunction{T}
    lower_bound::T
    pieces::Vector{AffineFunction{T}}
end
function function_value(f::PolyhedralFunction, x)
    return mapreduce(piece -> function_value(piece, x), max, f.pieces,
                     init = f.lower_bound)
end

function Base.:+(c::ConstantFunction, p::PolyhedralFunction)
    return PolyhedralFunction(c.value + p.lower_bound, p.pieces)
end

Base.:+(::ZeroFunction, f::Union{ConstantFunction, PolyhedralFunction}) = f

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
