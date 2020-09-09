struct HyperRectangle{VT<:AbstractArray}
    lb::VT
    ub::VT
end

# Changed all(f(x[i]) for i in eachindex(x)) to all(f.(x))
# See test_performances.
function Base.in(x::VT, rect::HyperRectangle{VT}) where VT
    return all(rect.lb .<= x .<= rect.ub)
end

function Base.isempty(rect::HyperRectangle)
    return any(rect.lb .> rect.ub)
end

function Base.intersect(a::S, b::S) where {S<:HyperRectangle}
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end

function Base.issubset(a::S, b::S) where {S<:HyperRectangle}
    return isempty(a) || (all(a.lb .>= b.lb) && all(a.ub .<= b.ub))
end

_ranges(rectI::HyperRectangle{SVector{N,Int}}) where N =
    ntuple(i -> UnitRange(rectI.lb[i], rectI.ub[i]), Val(N))

Base.eachindex(rect::HyperRectangle{<:SVector{N}}) where N = 1:N
