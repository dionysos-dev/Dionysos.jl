struct Ellipsoid{T <: Real, MT <: AbstractMatrix{T}, VT <: AbstractVector{T}}
    P::MT
    c::VT
    function Ellipsoid(
        P::MT,
        c::VT,
    ) where {T <: Real, MT <: AbstractMatrix{T}, VT <: AbstractVector{T}}
        M = (P + P') ./ 2
        if !isposdef(M) # see if delete this test
            error("matrix must be positive definite")
        end
        return new{T, MT, VT}(M, c)
    end
end

function get_center(elli::Ellipsoid)
    return elli.c
end

function get_shape(elli::Ellipsoid)
    return elli.P
end

function get_dims(elli::Ellipsoid)
    return length(elli.c)
end

function get_root(elli::Ellipsoid)
    return sqrt(elli.P)
end

function centerDistance(elli1::Ellipsoid, elli2::Ellipsoid)
    return norm(get_center(elli1) - get_center(elli2))
end

function pointCenterDistance(elli::Ellipsoid, x)
    return norm(get_center(elli) - x)
end

function get_volume(elli::Ellipsoid)
    N = size(elli.P, 1)
    return pi^(N / 2) / (gamma(N / 2 + 1)) * det(elli.P)^(-1 / 2)
end

function Base.:*(elli::Ellipsoid, r::Real)
    return Ellipsoid(elli.P * (1 / r), elli.c)
end

function Base.:*(r::Real, elli::Ellipsoid)
    return elli * r
end

function Base.:/(elli::Ellipsoid, r::Real)
    return elli * (1 / r)
end

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P * (1 / α), elli.c * α)
end

function expand(elli::Ellipsoid, α)
    return Ellipsoid(elli.P * (1 / α), elli.c)
end

function transform(elli::Ellipsoid, A, b)
    return Ellipsoid(A' \ elli.P / A, A * elli.c + b)
end

# return the ellidpoid f(Ε) where E = {x : (x-c)'P(x-c) <= 1} and f(x) = Ax+B
# with A invertible
function affine_transformation(elli::Ellipsoid, A, b)
    return Ellipsoid(A' \ elli.P / A, A * elli.c + b)
end

# get the point along the i largest axis
function get_axis_points(elli::Ellipsoid, i)
    specDecomp = eigen(elli.P)
    vals = specDecomp.values
    vectors = specDecomp.vectors

    # Trouver l'indice de la plus grande valeur propre
    sorted_indices = sortperm(vals)
    index = sorted_indices[i]
    # Extraire le vecteur propre correspondant
    vp = vectors[:, index]
    # Calculer la longueur de l'axe
    l = 1 / sqrt(vals[index])
    # Calculer les coordonnées des deux points
    p1 = elli.c - l * vp
    p2 = elli.c + l * vp
    return p1, p2
end

# get the point aling the axis, from the largest
function get_all_axis_points(elli::Ellipsoid)
    specDecomp = eigen(elli.P)
    vals = specDecomp.values
    vectors = specDecomp.vectors
    # Trouver l'indice de la plus grande valeur propre
    sorted_indices = sortperm(vals)
    L = []
    for i in 1:length(vals)
        index = sorted_indices[i]
        # Extraire le vecteur propre correspondant
        vp = vectors[:, index]
        # Calculer la longueur de l'axe
        l = 1 / sqrt(vals[index])
        # Calculer les coordonnées des deux points
        p1 = elli.c - l * vp
        p2 = elli.c + l * vp
        push!(L, (p1, p2))
    end
    return L
end

# get the radius of the largest ball that is inscribed in an ellipsoid of the i largest
# function with argument 1 returns the length of the longest semi-axis of the ellipsoid. 
function get_length_semiaxis_sorted(elli::Ellipsoid, i)
    specDecomp = eigen(elli.P)
    vals = specDecomp.values
    sorted_vals = sort(vals)
    return 1 / sqrt(sorted_vals[i])
end

function get_length_semiaxis(elli::Ellipsoid, i)
    specDecomp = eigen(elli.P)
    vals = specDecomp.values
    return 1 / sqrt(vals[i])
end

function get_length_semiaxis(elli::Ellipsoid)
    specDecomp = eigen(elli.P)
    vals = specDecomp.values
    return 1 ./ sqrt.(vals)
end

function get_inscribed_ball(elli::Ellipsoid)
    r = get_length_semiaxis_sorted(elli, 1)
    I_elli = Matrix{Float64}(I, size(elli.P)...)
    return Ellipsoid((1 / (r * r)) * I_elli, elli.c)
end

@recipe function f(e::Ellipsoid; axis_plot = false, color1 = :black, color2 = :black)
    if axis_plot
        @series begin
            color := color1
            p1, p2 = get_axis_points(e, 1)
            return DrawSegment(p1, p2)
        end
        color := color2
        p1, p2 = get_axis_points(e, 2)
        return DrawSegment(p1, p2)
    else
        opacity := 1.0
        label := ""
        lw := 1
        lc := :black
        Pvar = get_shape(e)
        Qvar = inv(Pvar)
        Qvar = (Qvar + Qvar') ./ 2
        return LazySets.Ellipsoid(collect(get_center(e)), Qvar)
    end
end

# get the farthest point of the ellipsoid in direction d
function get_farthest_point(elli::Ellipsoid, d)
    d = d / norm(d)
    Q = inv(elli.P)
    a = Q * d
    return a / sqrt(d'a)
end

"""
    get_min_bounding_box(elli, optimizer) 

Finds the minimum bounding box containing the ellipsoid {(x-c)'P(x-c) < 1}. 
"""
function get_min_bounding_box(elli::Ellipsoid; optimizer = nothing)
    P = elli.P
    n = size(P, 1)
    R = zeros(n)

    if optimizer !== nothing
        model = Model(optimizer)
        @variable(model, x[i = 1:n])
        @constraint(model, x'P * x <= 1)
        for i in 1:n
            new_model, reference_map = copy_model(model)
            set_optimizer(new_model, optimizer)
            @objective(new_model, Max, reference_map[x[i]])
            optimize!(new_model)
            R[i] = abs(value(reference_map[x[i]]))
        end
    else
        for i in 1:n
            ei = zeros(n)
            ei[i] = 1
            R[i] = get_farthest_point(elli, ei)[i]
        end
    end
    box = IntervalBox(elli.c .- R, elli.c .+ R)
    return box
end

function sample(elli::Ellipsoid; N = 500)
    box = get_min_bounding_box(elli)
    points = [sample(box) for i in 1:N]
    filter!(x -> x ∈ elli, points)
    return points
end

include("ellipsoid_inclusion.jl")
include("ellipsoid_intersection.jl")

# compress E1 if E1∩E2≠∅
# return nothing if impossible
function compress_if_intersection(E1::Ellipsoid, E2::Ellipsoid)
    if is_intersected(E1, E2)
        return scale_for_noninclusion_contact_point(E1, E2)
    else
        return E1
    end
end
