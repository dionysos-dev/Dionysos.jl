using LazySets

"""
    HyperRectangle{VT} <: LazySet

Defines a hyper-rectangle using `lb` (lower bound) and `ub` (upper bound).
"""
struct HyperRectangle{VT}
    lb::VT # Lower bounds of the rectangle
    ub::VT # Upper bounds of the rectangle
end

Base.in(x, rect::HyperRectangle) = all(rect.lb .<= x .<= rect.ub)
Base.in(rect1::HyperRectangle, rect2::HyperRectangle) =
    all(rect1.lb .>= rect2.lb) && all(rect1.ub .<= rect2.ub)
Base.isequal(rect1::HyperRectangle, rect2::HyperRectangle) =
    all(rect1.lb .== rect2.lb) && all(rect1.ub .== rect2.ub)
Base.:(==)(rect1::HyperRectangle, rect2::HyperRectangle) = isequal(rect1, rect2)
Base.isempty(rect::HyperRectangle) = any(rect.lb .> rect.ub)
is_intersection(a::HyperRectangle, b::HyperRectangle) = !Base.isempty(Base.intersect(a, b))
Base.intersect(a::HyperRectangle, b::HyperRectangle) =
    HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
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

function get_vertices(rect::HyperRectangle)
    n = length(rect.lb)
    vertices = zeros(n, 2^n)
    for i in 0:(2^n - 1)
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

function rectangle(c, r)
    return Shape(
        c[1] .- r[1] .+ [0, 2 * r[1], 2 * r[1], 0],
        c[2] .- r[2] .+ [0, 0, 2 * r[2], 2 * r[2]],
    )
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

"""
    DeformedRectangleDraw(rec, N, f)

A helper struct for drawing a deformed rectangle.
"""
struct DeformedRectangleDraw
    rect::HyperRectangle
    f::Function
    N::Integer # used for precision
    shape::Any
end

function DeformedRectangleDraw(rect::HyperRectangle, f::Function; N = 10, dims = [1, 2])
    lb = rect.lb[dims]
    ub = rect.ub[dims]
    points = SVector[]
    for x in LinRange(lb[1], ub[1], N)
        push!(points, f(SVector(x, lb[2])))
    end
    for x in LinRange(lb[2], ub[2], N)
        push!(points, f(SVector(ub[1], x)))
    end
    for x in LinRange(ub[1], lb[1], N)
        push!(points, f(SVector(x, ub[2])))
    end
    for x in LinRange(ub[2], lb[2], N)
        push!(points, f(SVector(lb[1], x)))
    end
    unique!(points)
    x = [point[1] for point in points]
    y = [point[2] for point in points]

    return DeformedRectangleDraw(rect, f, N, Shape(x, y))
end

@recipe function f(deformed_rect::DeformedRectangleDraw; dims = [1, 2])
    @series begin
        return deformed_rect.shape
    end
end
