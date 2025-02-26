using LazySets

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
    return isequal(rect1, rect2)
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
        return prod(rect.ub - rect.lb)
    end
end

function is_intersection(a::HyperRectangle, b::HyperRectangle)
    return !Base.isempty(Base.intersect(a, b))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
end

function get_center(rect::HyperRectangle)
    return (rect.lb + rect.ub) / 2
end

function get_h(rect::HyperRectangle)
    return rect.ub - rect.lb
end

function get_r(rect::HyperRectangle)
    return get_h(rect) ./ 2.0
end

function get_dims(rect::HyperRectangle)
    return length(rect.lb)
end

function scale(rect::HyperRectangle, α)
    return HyperRectangle(rect.lb * α, rect.ub * α)
end

function get_volume(rect::HyperRectangle)
    return prod(rect.ub - rect.lb)
end

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
    n = get_dims(rect)
    sample_point = similar(rect.lb)
    for i in 1:n
        sample_point[i] = rand() * (rect.ub[i] - rect.lb[i]) + rect.lb[i]
    end
    return sample_point
end

function affine_transformation(rect::HyperRectangle, A, b)
    X = LazySets.Hyperrectangle(Vector(get_center(rect)), Vector(get_r(rect)))
    return LazySets.AffineMap(Matrix(A), X, Vector(b))
end

function rectangle(c, r)
    return Shape(
        c[1] .- r[1] .+ [0, 2 * r[1], 2 * r[1], 0],
        c[2] .- r[2] .+ [0, 0, 2 * r[2], 2 * r[2]],
    )
end

RecipesBase.@recipe function f(rect::HyperRectangle; dims = [1, 2])
    center = get_center(rect)
    h = get_h(rect)
    return rectangle(center[dims], h[dims] ./ 2)
end
