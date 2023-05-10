using SpecialFunctions
using Plots
using LazySets
using LinearAlgebra
using Polyhedra
using IntervalArithmetic
using JuMP

# model (x-c)' U' U (x-c) ≤ 1 <=> ||U(x*c)||≤1
# note that the differnece with Ellipdoid structure is that the matrix U'u
# could not be invertible
struct DegenerateEllipsoid{T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    U::MT
    c::VT
    P::MT
    function DegenerateEllipsoid(U::MT,c::VT) where {T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
        return new{T,MT,VT}(U,c,U'*U)
    end
end

function get_center(elli::DegenerateEllipsoid)
    return elli.c
end

function get_shape(elli::DegenerateEllipsoid)
    return elli.P
end

function get_dims(elli::DegenerateEllipsoid)
    return length(elli.c)
end

function affine_transformation(elli::DegenerateEllipsoid, A, b)
    P = get_shape(elli)
    return Ellipsoid(A'\P/A, A*elli.c+b) 
end

function get_root(elli::DegenerateEllipsoid)
    return elli.U
end

function is_degenerate(elli::DegenerateEllipsoid)
    P = get_shape(elli)
    return !isposdef((P+P')./2)
end

# If U'U is not invertible, the ellipsoid defined by {x: (x-c)'U'U(x-c) <= 1} is a 
# degenerate ellipsoid and can be visualized as a line segment
function plotE!(elli::DegenerateEllipsoid; color=:blue, opacity=1.0, label="",lw=1,lc=:black)
    P = get_shape(elli)
    if is_degenerate(elli)
        U = get_root(elli)
        v = nullspace(U'*U)
        v = v ./ norm(v)
        c = get_center(elli)
        # Define two points on the line segment
        p1 = c - v
        p2 = c + v

        # Plot the line segment
        plot([p1[1], p2[1]], [p1[2], p2[2]], xlims=(minimum([p1[1], p2[1]])-0.1, maximum([p1[1], p2[1]])+0.1), ylims=(minimum([p1[2], p2[2]])-0.1, maximum([p1[2], p2[2]])+0.1), aspect_ratio=:equal)
    else
        Q = inv(P)
        Q = (Q+Q')./2
        E = LazySets.Ellipsoid(collect(get_center(elli)), Q) #not optimal, require to inverse
        Plots.plot!(E; color=color, opacity=opacity, label=label,lw=lw,lc=lc)
    end
end