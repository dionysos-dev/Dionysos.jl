
_normalize_dims(dims) =
    dims isa Tuple{Int, Int} ? dims :
    dims isa AbstractVector{<:Integer} ? (Int(dims[1]), Int(dims[2])) :
    throw(ArgumentError("dims must be (i,j) or [i,j], got $(typeof(dims))"))

struct DrawPoint{T <: Real, VT <: AbstractVector{T}}
    p::VT
end
@recipe function f(p::DrawPoint; dims = [1, 2])
    i, j = _normalize_dims(dims)
    color --> :black
    marker --> :circle
    linetype --> :scatter
    label := ""
    return [p.p[i]], [p.p[j]]
end

struct DrawArrow{T <: Real, VT <: AbstractVector{T}}
    p1::DrawPoint{T, VT}
    p2::DrawPoint{T, VT}
    function DrawArrow(
        p1::DrawPoint{T, VT},
        p2::DrawPoint{T, VT},
    ) where {T <: Real, VT <: AbstractVector{T}}
        return new{T, VT}(p1, p2)
    end
    function DrawArrow(p1::VT, p2::VT) where {T <: Real, VT <: AbstractVector{T}}
        return new{T, VT}(DrawPoint(p1), DrawPoint(p2))
    end
end
@recipe function f(a::DrawArrow; dims = [1, 2])
    i, j = _normalize_dims(dims)
    marker --> :circle
    markeralpha --> 0.0
    color --> :black
    arrow --> (:closed, 2.0)
    label := ""
    return [a.p1.p[i], a.p2.p[i]], [a.p1.p[j], a.p2.p[j]]
end

struct DrawSegment{T <: Real, VT <: AbstractVector{T}}
    p1::DrawPoint{T, VT}
    p2::DrawPoint{T, VT}
    function DrawSegment(
        p1::DrawPoint{T, VT},
        p2::DrawPoint{T, VT},
    ) where {T <: Real, VT <: AbstractVector{T}}
        return new{T, VT}(p1, p2)
    end
    function DrawSegment(p1::VT, p2::VT) where {T <: Real, VT <: AbstractVector{T}}
        return new{T, VT}(DrawPoint(p1), DrawPoint(p2))
    end
end
@recipe function f(s::DrawSegment, dims = [1, 2])
    i, j = _normalize_dims(dims)
    linestyle --> :dash
    color --> :black
    label := ""
    return [s.p1.p[i], s.p2.p[i]], [s.p1.p[j], s.p2.p[j]]
end

struct DrawTrajectory{
    T <: Real,
    VT <: AbstractVector{T},
    VP <: AbstractVector{DrawPoint{T, VT}},
}
    vp::VP
    function DrawTrajectory(
        vp::VP,
    ) where {T <: Real, VT <: AbstractVector{T}, VP <: AbstractVector{DrawPoint{T, VT}}}
        return new{T, VT, VP}(vp)
    end
    function DrawTrajectory(
        trajCoord::VVT,
    ) where {T <: Real, VT <: AbstractVector{T}, VVT <: AbstractVector{VT}}
        vec = [DrawPoint(coord) for coord in trajCoord]
        VP = typeof(vec)
        return new{T, VT, VP}(vec)
    end
end

@recipe function f(t::DrawTrajectory; dims = [1, 2])
    for i in 1:(length(t.vp) - 1)
        @series begin
            dims := dims
            t.vp[i + 1]
        end
        @series begin
            dims := dims
            DrawArrow(t.vp[i], t.vp[i + 1])
        end
    end
    @series begin
        dims := dims
        t.vp[1]
    end
end
