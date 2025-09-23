
struct DrawPoint{T <: Real, VT <: AbstractVector{T}}
    p::VT
end
@recipe function f(p::DrawPoint)
    color --> :black
    marker --> :circle
    linetype --> :scatter
    legend --> false
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
    arrow --> (:closed, 1.0)
    legend --> false
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
    legend --> false
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
    ##solve finding center (other solvers? https://jump.dev/JuMP.jl/dev/installation/#Supported-solvers)
    @suppress begin
        solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)
        if t !== nothing
            c, r = hchebyshevcenter(hrep(po), solver; verbose = 0)
            annotate!(fig, c[1], c[2], t)
        end
    end
end

struct DrawCostTrajectory{
    T <: Real,
    VT <: AbstractVector{T},
    VP <: AbstractVector{DrawPoint{T, VT}},
    VC <: AbstractVector{Any},
}
    vp::VP
    costs::VC
    function DrawCostTrajectory(
        vp::VP, costs::VC
    ) where {T <: Real, VT <: AbstractVector{T}, VP <: AbstractVector{DrawPoint{T, VT}}, VC <: AbstractVector{Any}}
        return new{T, VT, VP, VC}(vp, costs)
    end
    function DrawCostTrajectory(
        trajcoord::VVT, costs::VC
    ) where {T <: Real, VT <: AbstractVector{T}, VVT <: AbstractVector{VT}, VC <: AbstractVector{Any}}
        vec = [DrawPoint(coord) for coord in trajcoord]
        VP = typeof(vec)
        return new{T, VT, VP, VC}(vec, costs)
    end 
end

@recipe function f(traj::DrawCostTrajectory)
    min_cost, max_cost = extrema(traj.costs) 
    println(sort(unique(traj.costs)))
    colormap = Colors.colormap("Oranges") #reverse(cgrad(:blues))
    println("$min_cost, $max_cost")

    for i in 1:(length(traj.vp) - 1)
        normalized_cost = (traj.costs[i] - min_cost) / (max_cost - min_cost + eps())  
        normalized_cost = round(Int, normalized_cost * 100)
        color = colormap[min(100, max(1, normalized_cost))]
        @series begin
            seriescolor --> color
            DrawArrow(traj.vp[i], traj.vp[i + 1])
        end
    end
    return traj.vp[1] 
end