export ZeroFunction, ConstantFunction, QuadraticControlFunction, PolyhedralFunction
export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export OptimalControlProblem
export AffineFunction, PolyhedralFunction
export last_mode, function_value

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
function DiscreteTrajectory{TT}(q_0) where TT
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

"""
    ContinuousTrajectory{T, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}

`x` is a sequence of points in the state space and `u` is a sequence of points
in the input space
"""
struct ContinuousTrajectory{T, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}
    x::Vector{XVT}
    u::Vector{UVT}
end

struct ContinuousTrajectoryAttribute <: MOI.AbstractModelAttribute
end

struct HybridTrajectory{T, TT, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T,XVT,UVT}
end

struct ZeroFunction end

struct ConstantFunction{T}
    value::T
end
function_value(f::ConstantFunction, x) = f.value
function Base.:+(f::ConstantFunction, g::ConstantFunction)
    return ConstantFunction(f.value + g.value)
end

struct QuadraticControlFunction{T, MT<:AbstractMatrix{T}}
    Q::MT
end

struct AffineFunction{T}
    a::Vector{T}
    β::T
end


function function_value(f::AffineFunction, x)
    return sum(f.a .* x) + f.β
end
function Base.isapprox(f::AffineFunction, g::AffineFunction; kws...)
    return isapprox(f.a, g.a; kws...) && isapprox(f.β, g.β; kws...)
end

struct PolyhedralFunction{T}
    lower_bound::T
    pieces::Vector{AffineFunction{T}}
    domain::Polyhedra.Intersection{T,Vector{T},Int}
end
_inf(T::Type{<:AbstractFloat}) = typemax(T)
_inf(T::Type) = error("No infinite value for type $T")
function function_value(f::PolyhedralFunction{T}, x) where T
    if !(x in f.domain)
        return _inf(T)
    end
    return mapreduce(piece -> function_value(piece, x), max,
        f.pieces, init = f.lower_bound)
end

function Base.:+(c::ConstantFunction, p::PolyhedralFunction)
    return PolyhedralFunction(c.value + p.lower_bound, p.pieces, p.domain)
end

Base.:+(::ZeroFunction, f::Union{ConstantFunction,PolyhedralFunction}) = f

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
