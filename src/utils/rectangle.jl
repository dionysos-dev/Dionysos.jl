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

function isintersect(a::HyperRectangle, b::HyperRectangle)
    return !Base.isempty(Base.intersect(a, b))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
end
