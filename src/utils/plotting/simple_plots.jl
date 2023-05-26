using Plots

struct DrawPoint{T<:Real, VT<:AbstractVector{T}}
    p::VT
end
@recipe function f(p::DrawPoint)
    color       --> :black
    marker      --> :circle
    linetype    --> :scatter
    legend      --> false
    [p.p[1]], [p.p[2]]
end

struct DrawArrow{T<:Real, VT<:AbstractVector{T}}
    p1::DrawPoint{T, VT}
    p2::DrawPoint{T, VT}
    function DrawArrow(p1::DrawPoint{T, VT}, p2::DrawPoint{T, VT}) where {T<:Real, VT<:AbstractVector{T}}
        new{T, VT}(p1, p2)
    end
    function DrawArrow(p1::VT, p2::VT) where {T<:Real, VT<:AbstractVector{T}}
        new{T, VT}(DrawPoint(p1), DrawPoint(p2))
    end
end
@recipe function f(a::DrawArrow) 
    marker      --> :circle
    markeralpha --> 0.
    color       --> :black
    arrow       --> (:closed, 2.)
    legend      --> false
    [a.p1.p[1], a.p2.p[1]], [a.p1.p[2], a.p2.p[2]] 
end

struct DrawSegment{T<:Real, VT<:AbstractVector{T}}
    p1::DrawPoint{T, VT}
    p2::DrawPoint{T, VT}
    function DrawSegment(p1::DrawPoint{T, VT}, p2::DrawPoint{T, VT}) where {T<:Real, VT<:AbstractVector{T}}
        new{T, VT}(p1, p2)
    end
    function DrawSegment(p1::VT, p2::VT) where {T<:Real, VT<:AbstractVector{T}}
        new{T, VT}(DrawPoint(p1), DrawPoint(p2)) 
    end
end
@recipe function f(s::DrawSegment)
    linestyle   --> :dash
    color       --> :black
    legend      --> false
    [s.p1.p[1], s.p2.p[1]], [s.p1.p[2], s.p2.p[2]]
end

struct DrawTrajectory{T<:Real, VT<:AbstractVector{T}, VP<:AbstractVector{DrawPoint{T, VT}}}
    vp::VP
    function DrawTrajectory(vp::VP) where {T<:Real, VT<:AbstractVector{T}, VP<:AbstractVector{DrawPoint{T, VT}}}
        new{T, VT, VP}(vp)
    end
    function DrawTrajectory(trajCoord::VVT) where {T<:Real, VT<:AbstractVector{T}, VVT<:AbstractVector{VT}}
        vec = [DrawPoint(coord) for coord in trajCoord]
        VP = typeof(vec)
        new{T, VT, VP}(vec)
    end
end

@recipe function f(t::DrawTrajectory)
    for i = 1:length(t.vp) - 1
        @series begin t.vp[i + 1] end
        @series begin DrawArrow(t.vp[i], t.vp[i + 1]) end
    end
    t.vp[1]
end
