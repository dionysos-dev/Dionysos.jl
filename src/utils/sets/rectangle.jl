import Base

"""
    HyperRectangle{N,T}

Axis-aligned hyper-rectangle with lower bound `lb` and upper bound `ub`.
"""
struct HyperRectangle{N,T} <: AbstractSetNode{N,T}
    lb::SVector{N,T}
    ub::SVector{N,T}
end

# From tuples
HyperRectangle(lb::NTuple{N,Ti}, ub::NTuple{N,Ti}) where {N,Ti} = HyperRectangle(SVector{N,Ti}(lb), SVector{N,Ti}(ub))

# From vectors (runtime dimension)
HyperRectangle(lb::AbstractVector{Ti}, ub::AbstractVector{Ti}) where {Ti} = HyperRectangle(SVector{length(lb),Ti}(lb), SVector{length(ub),Ti}(ub))

Base.in(x, rect::HyperRectangle) = all(rect.lb .<= x .<= rect.ub)
Base.in(rect1::HyperRectangle, rect2::HyperRectangle) =
    all(rect1.lb .>= rect2.lb) && all(rect1.ub .<= rect2.ub)
Base.isequal(rect1::HyperRectangle, rect2::HyperRectangle) =
    all(rect1.lb .== rect2.lb) && all(rect1.ub .== rect2.ub)
Base.:(==)(rect1::HyperRectangle, rect2::HyperRectangle) = isequal(rect1, rect2)
Base.isempty(rect::HyperRectangle) = any(rect.lb .> rect.ub)
is_intersection(a::HyperRectangle, b::HyperRectangle) = !Base.isempty(Base.intersect(a, b))
import Base: intersect
function intersect(a::HyperRectangle{N,T}, b::HyperRectangle{N,T}) where {N,T}
    HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end
Base.issubset(a::HyperRectangle, b::HyperRectangle) =
    all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
get_center(rect::HyperRectangle) = (rect.lb + rect.ub) / 2
get_h(rect::HyperRectangle) = rect.ub - rect.lb
get_r(rect::HyperRectangle) = get_h(rect) ./ 2.0
get_dims(rect::HyperRectangle) = length(rect.lb)
volume(rect::HyperRectangle) = Base.isempty(rect) ? 0.0 : prod(rect.ub - rect.lb)
scale(rect::HyperRectangle, α) = HyperRectangle(rect.lb * α, rect.ub * α)
to_LazySets(rect::HyperRectangle) =
    LazySets.Hyperrectangle(Vector(get_center(rect)), Vector(get_r(rect)))
affine_transformation(rect::HyperRectangle, A, b) =
    LazySets.AffineMap(Matrix(A), to_LazySets(rect), Vector(b))
get_volume(rect::HyperRectangle) = prod(rect.ub .- rect.lb)

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
    rect::HyperRectangle{N,T},
    lb::NTuple{N,T},
    ub::NTuple{N,T},
    periodic_dims::SVector{P,Int},
    periods::SVector{P,T},
    start::SVector{P,T},
    i::Int,
) where {N,T,P}
    if i > P
        push!(out, HyperRectangle(SVector{N,T}(lb), SVector{N,T}(ub)))
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
    set_in_period(rect, periodic_dims, periods, start) -> LazySetUnion{N,T}

Split `rect` along periodic boundaries and return a `LazySetUnion` of wrapped rectangles.
"""
function set_in_period(
    rect::HyperRectangle{N,T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
) where {N,T,P}
    L = HyperRectangle{N,T}[]
    _recursive_period_split!(L, rect, Tuple(rect.lb), Tuple(rect.ub),
                            periodic_dims, periods, start, 1)

    # you need a constructor; otherwise build explicitly:
    U = LazySetUnion{N,T}()
    Base.append!(U.sets, L)
    return U
end

"""
    DeformedRectangle(rect, f, N, shape)

A deformed rectangle.
"""
struct DeformedRectangle{N,T} <: AbstractSetNode{N,T}
    rect::HyperRectangle{N,T}
    f::Function
end

function SampleBoundaryDeformedRectangle(drect::DeformedRectangle{N,T}; K=50, dims=(1,2)) where {N,T}
    rect = drect.rect
    f = drect.f
    d1, d2 = dims
    lb = rect.lb[[d1,d2]]
    ub = rect.ub[[d1,d2]]

    points = SVector{2,T}[]  # assumes dims=(1,2)

    for x in LinRange(lb[1], ub[1], K)
        push!(points, f(SVector(x, lb[2])))
    end
    for y in LinRange(lb[2], ub[2], K)
        push!(points, f(SVector(ub[1], y)))
    end
    for x in LinRange(ub[1], lb[1], K)
        push!(points, f(SVector(x, ub[2])))
    end
    for y in LinRange(ub[2], lb[2], K)
        push!(points, f(SVector(lb[1], y)))
    end

    return points
end

@recipe function f(drect::DeformedRectangle; dims=(1,2), K=50)
    pts = SampleBoundaryDeformedRectangle(drect; K=K, dims=dims)
    x = getindex.(pts, 1)
    y = getindex.(pts, 2)
    @series begin
        seriestype := :polygon
        (x, y)
    end
end