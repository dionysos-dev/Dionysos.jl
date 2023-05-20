struct HyperRectangle{VT}
    lb::VT
    ub::VT
end

# Changed all(f(x[i]) for i in eachindex(x)) to all(f.(x))
# See test_performances.
function Base.in(x, rect::HyperRectangle)
    return all(rect.lb .<= x .<= rect.ub)
end

function Base.in(rect1::HyperRectangle, rect2::HyperRectangle)
    return all(rect1.lb .>= rect2.lb) && all(rect1.ub .<= rect2.ub)
end

function Base.isequal(rect1::HyperRectangle, rect2::HyperRectangle)
    return all(rect1.lb .== rect2.lb) && all(rect1.ub .== rect2.ub)
end

function Base.:(==)(rect1::HyperRectangle, rect2::HyperRectangle)
    return isequal(rect1,rect2)    
end
function Base.isempty(rect::HyperRectangle)
    return any(rect.lb .> rect.ub)
end

function Base.intersect(a::HyperRectangle, b::HyperRectangle)
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end

function volume(rect::HyperRectangle)
    if Base.isempty(rect)
        return 0.0
    else
        return prod(rect.ub-rect.lb)
    end
end

function is_intersection(a::HyperRectangle, b::HyperRectangle)
    return !Base.isempty(Base.intersect(a, b))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
end

function get_center(rect::HyperRectangle)
    return (rect.lb+rect.ub)/2
end

function get_h(rect::HyperRectangle)
    return rect.ub-rect.lb
end

function get_dims(rect::HyperRectangle)
    return length(rect.lb)
end

function scale(rect::HyperRectangle, α)
    return HyperRectangle(rect.lb*α,rect.ub*α)
end

function rectangle(c,r)
    Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
end

function Plots.plot!(rect::HyperRectangle; dims=[1,2], color=:yellow, opacity=1.0)
    center = get_center(rect)
    h = get_h(rect)
    plot!(rectangle(center[dims],h[dims]./2), opacity=opacity, color=color)
end