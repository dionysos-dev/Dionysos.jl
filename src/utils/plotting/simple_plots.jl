
struct DrawPoint{T <: Real, VT <: AbstractVector{T}}
    p::VT
end
@recipe function f(p::DrawPoint)
    color --> :black
    marker --> :circle
    linetype --> :scatter
    label := ""
    return [p.p[1]], [p.p[2]]
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
@recipe function f(a::DrawArrow)
    marker --> :circle
    markeralpha --> 0.0
    color --> :black
    arrow --> (:closed, 2.0)
    label := ""
    return [a.p1.p[1], a.p2.p[1]], [a.p1.p[2], a.p2.p[2]]
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
@recipe function f(s::DrawSegment)
    linestyle --> :dash
    color --> :black
    label := ""
    return [s.p1.p[1], s.p2.p[1]], [s.p1.p[2], s.p2.p[2]]
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

@recipe function f(t::DrawTrajectory)
    for i in 1:(length(t.vp) - 1)
        @series begin
            t.vp[i + 1]
        end
        @series begin
            DrawArrow(t.vp[i], t.vp[i + 1])
        end
    end
    return t.vp[1]
end

# Auxiliary function for annotation
function text_in_set_plot!(fig, po::Polyhedra.Rep, t)
    ##solve finding center (HiGHS is currently the best open source LP solver)
    solver = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    if t !== nothing
        c, r = Polyhedra.hchebyshevcenter(Polyhedra.hrep(po), solver; verbose = 0)
        annotate!(fig, c[1], c[2], t)
    end
end
