struct HyperRectangle{VT}
    lb::VT
    ub::VT
end

# Changed all(f(x[i]) for i in eachindex(x)) to all(f.(x))
# See test_performances.
function Base.in(x, rect::HyperRectangle)
    return all(rect.lb .<= x .<= rect.ub)
end

function Base.isempty(rect::HyperRectangle)
    return any(rect.lb .> rect.ub)
end

function Base.intersect(a::HyperRectangle, b::HyperRectangle)
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
end
