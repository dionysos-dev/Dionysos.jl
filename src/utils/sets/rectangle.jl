import Base

"""
    HyperRectangle{N,T}

Axis-aligned hyper-rectangle with lower bound `lb` and upper bound `ub`.
"""
struct HyperRectangle{N, T} <: AbstractSetNode{N, T}
    lb::SVector{N, T}
    ub::SVector{N, T}
end

# From tuples
HyperRectangle(lb::T, ub::T) where {T <: Tuple} = HyperRectangle(SVector(lb), SVector(ub))

# From vectors (runtime dimension)
function HyperRectangle(lb::AbstractVector{Ti}, ub::AbstractVector{Ti}) where {Ti}
    n = length(lb)
    # Guard: without this, a length mismatch keeps re-matching this same
    # AbstractVector method (the `SVector{length,·}` results are still vectors of
    # unequal length) and recurses until a StackOverflowError.
    n == length(ub) || throw(
        DimensionMismatch(
            "lb and ub must have equal length, got $(length(lb)) and $(length(ub))",
        ),
    )
    return HyperRectangle(SVector{n, Ti}(lb), SVector{n, Ti}(ub))
end

Base.in(x::AbstractVector, rect::HyperRectangle) = all(rect.lb .<= x .<= rect.ub)
Base.isequal(rect1::HyperRectangle, rect2::HyperRectangle) =
    all(rect1.lb .== rect2.lb) && all(rect1.ub .== rect2.ub)
Base.:(==)(rect1::HyperRectangle, rect2::HyperRectangle) = isequal(rect1, rect2)
Base.isempty(rect::HyperRectangle) = any(rect.lb .> rect.ub)
Base.isdisjoint(a::HyperRectangle, b::HyperRectangle) = Base.isempty(Base.intersect(a, b))
import Base: intersect
function intersect(a::HyperRectangle{N, T}, b::HyperRectangle{N, T}) where {N, T}
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end
Base.issubset(a::HyperRectangle, b::HyperRectangle) =
    all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
get_center(rect::HyperRectangle) = (rect.lb + rect.ub) / 2
get_h(rect::HyperRectangle) = rect.ub - rect.lb
get_r(rect::HyperRectangle) = get_h(rect) ./ 2.0
get_dim(rect::HyperRectangle) = length(rect.lb)

LazySets.center(rect::HyperRectangle) = get_center(rect)
LazySets.low(rect::HyperRectangle) = rect.lb
LazySets.high(rect::HyperRectangle) = rect.ub
LazySets.low(rect::HyperRectangle, i::Int) = rect.lb[i]
LazySets.high(rect::HyperRectangle, i::Int) = rect.ub[i]
LazySets.ρ(d::AbstractVector, rect::HyperRectangle) =
    d ⋅ get_center(rect) + abs.(d) ⋅ get_r(rect)
LazySets.σ(d::AbstractVector, rect::HyperRectangle) =
    get_center(rect) .+ get_r(rect) .* sign.(d)
scale(rect::HyperRectangle, α) = HyperRectangle(rect.lb * α, rect.ub * α)
to_LazySets(rect::HyperRectangle) =
    LazySets.Hyperrectangle(Vector(get_center(rect)), Vector(get_r(rect)))
affine_transformation(rect::HyperRectangle, A, b) =
    LazySets.AffineMap(Matrix(A), to_LazySets(rect), Vector(b))
get_volume(rect::HyperRectangle) =
    Base.isempty(rect) ? zero(eltype(rect.lb)) : prod(rect.ub .- rect.lb)

_outer_box(X::HyperRectangle) = X

function enlarge(rect::HyperRectangle, δ)
    return HyperRectangle(rect.lb .- δ, rect.ub .+ δ)
end

function get_vertices(rect::HyperRectangle)
    n = length(rect.lb)
    vertices = zeros(n, 2^n)
    for i in 0:(2 ^ n - 1)
        vertex = zeros(n)
        for j in 1:n
            vertex[j] = i & (1 << (j - 1)) > 0 ? rect.ub[j] : rect.lb[j]
        end
        vertices[:, i + 1] = vertex
    end
    return vertices
end

function collect_vertices(rect::HyperRectangle)
    vertices_matrix = get_vertices(rect)
    vertices = [vertices_matrix[:, i] for i in 1:size(vertices_matrix, 2)]
    return vertices
end

function sample(rect::HyperRectangle)
    return SVector{length(rect.lb)}(rand() .* (rect.ub .- rect.lb) .+ rect.lb)
end
function samples(rect::HyperRectangle, N::Int)
    return [sample(rect) for _ in 1:N]
end

@recipe function f(rect::HyperRectangle; dims = [1, 2])
    center = get_center(rect)[dims]
    r = get_h(rect)[dims] ./ 2  # Half-width
    x1, y1 = center[1] - r[1], center[2] - r[2]
    x2, y2 = center[1] + r[1], center[2] + r[2]

    @series begin
        seriestype := :shape
        [x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1]
    end
end

# --------------------------
# Periodic splitting
# --------------------------

function one_direction(lb, ub, Tper, T0)
    if ub - lb >= Tper
        return [(T0, T0 + Tper)]
    else
        lbw = T0 + mod(lb - T0, Tper)
        ubw = T0 + mod(ub - T0, Tper)
        if lbw <= ubw
            return [(lbw, ubw)]
        else
            return [(T0, ubw), (lbw, T0 + Tper)]
        end
    end
end

function _recursive_period_split!(
    out,
    rect::HyperRectangle{N, T},
    lb::NTuple{N, T},
    ub::NTuple{N, T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    i::Int,
) where {N, T, P}
    if i > P
        push!(out, HyperRectangle(SVector{N, T}(lb), SVector{N, T}(ub)))
        return
    end
    dim = periodic_dims[i]
    intervals = one_direction(rect.lb[dim], rect.ub[dim], periods[i], start[i])
    for (lI, uI) in intervals
        l = ntuple(k -> k == dim ? lI : lb[k], Val(N))
        u = ntuple(k -> k == dim ? uI : ub[k], Val(N))
        _recursive_period_split!(out, rect, l, u, periodic_dims, periods, start, i + 1)
    end
end

"""
    set_in_period(rect, periodic_dims, periods, start) -> LazySets.UnionSetArray

Split `rect` along periodic boundaries and return the union of the wrapped rectangles.
"""
function set_in_period(
    rect::HyperRectangle{N, T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
) where {N, T, P}
    L = HyperRectangle{N, T}[]
    _recursive_period_split!(
        L,
        rect,
        Tuple(rect.lb),
        Tuple(rect.ub),
        periodic_dims,
        periods,
        start,
        1,
    )
    return set_union(L)
end
