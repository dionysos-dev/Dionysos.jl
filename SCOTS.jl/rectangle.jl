struct HyperRectangle{T, VT<:AbstractVector{T}}
    lb::VT
    ub::VT
end
function Base.in(x::AbstractVector, rect::HyperRectangle)
    return all(i -> rect.lb[i] .<= x[i] .<= rect.ub[i], eachindex(x))
end
function Base.isempty(rect::HyperRectangle)
    return any(i -> rect.lb[i] > rect.ub[i], eachindex(rect.lb))
end
function Base.intersect(a::HyperRectangle, b::HyperRectangle)
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end
function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(i -> a.lb[i] >= b.lb[i], eachindex(a.lb)) &&
        all(i -> a.ub[i] <= b.ub[i], eachindex(a.ub))
end
