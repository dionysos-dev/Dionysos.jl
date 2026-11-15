# model (x-c)' U' U (x-c) ≤ 1 <=> ||U(x*c)||≤1
# note that the differnece with Ellipdoid structure is that the matrix U'u
# could not be invertible
struct DegenerateEllipsoid{T <: Real, MT <: AbstractMatrix{T}, VT <: AbstractVector{T}}
    U::MT
    c::VT
    P::MT
    function DegenerateEllipsoid(
        U::MT,
        c::VT,
    ) where {T <: Real, MT <: AbstractMatrix{T}, VT <: AbstractVector{T}}
        return new{T, MT, VT}(U, c, U' * U)
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
    return Ellipsoid(A' \ P / A, A * elli.c + b)
end

function get_root(elli::DegenerateEllipsoid)
    return elli.U
end

function is_degenerate(elli::DegenerateEllipsoid)
    P = get_shape(elli)
    return !isposdef((P + P') ./ 2)
end

# If U'U is not invertible, the ellipsoid defined by {x: (x-c)'U'U(x-c) <= 1} is a 
# degenerate ellipsoid and can be visualized as a line segment
@recipe function f(e::DegenerateEllipsoid)
    if !is_degenerate(e)
        return Ellipsoid(e.P, e.c)
    end

    color := :blue
    opacity := 1.0
    label := ""
    lw := 1
    lc := :black

    U = get_root(e)
    v = nullspace(U' * U)
    v = v ./ norm(v)
    c = get_center(e)
    # Define two points on the line segment
    p1 = c - v
    p2 = c + v

    xlims := (minimum([p1[1], p2[1]]) - 0.1, maximum([p1[1], p2[1]]) + 0.1)
    ylims := (minimum([p1[2], p2[2]]) - 0.1, maximum([p1[2], p2[2]]) + 0.1)
    aspect_ratio := equal

    return DrawSegment(p1, p2)
end
