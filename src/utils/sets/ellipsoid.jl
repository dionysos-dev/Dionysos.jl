struct Ellipsoid{N, T, MT <: SMatrix{N, N, T}, VT <: SVector{N, T}} <: AbstractSetNode{N, T}
    P::MT
    c::VT
    function Ellipsoid(
        P::MT,
        c::VT,
    ) where {N, T, MT <: SMatrix{N, N, T}, VT <: SVector{N, T}}
        M = (P + P')/2
        isposdef(Matrix(M)) || error("matrix must be positive definite")
        return new{N, T, MT, VT}(M, c)
    end
end

function Ellipsoid(P::AbstractMatrix{T}, c::AbstractVector{T}) where {T <: Real}
    n, m = size(P)
    n == m || throw(ArgumentError("P must be square, got $(size(P))"))
    length(c) == n || throw(ArgumentError("c must have length $n, got $(length(c))"))

    PS = SMatrix{n, n, T}(P)
    cS = SVector{n, T}(c)
    return Ellipsoid(PS, cS)
end

get_center(e::Ellipsoid) = e.c
get_shape(elli::Ellipsoid) = elli.P

function get_dim(elli::Ellipsoid)
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

# Volume of {x : (x-c)'P(x-c) ≤ 1} = π^(N/2) / Γ(N/2 + 1) · det(P)^(-1/2)
# (the volume of the unit N-ball scaled by the ellipsoid's semi-axes).
function get_volume(elli::Ellipsoid)
    N = size(elli.P, 1)
    return pi^(N / 2) / SpecialFunctions.gamma(N / 2 + 1) * det(elli.P)^(-1 / 2)
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

_outer_box(elli::Ellipsoid) = get_min_bounding_box(elli)

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P * (1 / α), elli.c * α)
end

function get_sublevel_set(elli::Ellipsoid, α)
    return Ellipsoid(elli.P * (1 / α), elli.c)
end

# return the ellipsoid f(Ε) where E = {x : (x-c)'P(x-c) <= 1} and f(x) = Ax+b
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

@recipe function f(
    e::Ellipsoid{N, T};
    axis_plot = false,
    color1 = :black,
    color2 = :black,
) where {N, T}
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
        return LazySets.Ellipsoid(Vector(get_center(e)), Matrix(Qvar))
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
    get_min_bounding_box(elli)

Finds the minimum bounding box containing the ellipsoid {(x-c)'P(x-c) < 1}.
"""
function get_min_bounding_box(elli::Ellipsoid)
    P = elli.P
    n = size(P, 1)
    R = zeros(n)
    for i in 1:n
        ei = zeros(n)
        ei[i] = 1
        R[i] = get_farthest_point(elli, ei)[i]
    end
    return HyperRectangle(elli.c .- R, elli.c .+ R)
end

function sample(elli::Ellipsoid; N = 500)
    rec = get_min_bounding_box(elli)
    points = [sample(rec) for i in 1:N]
    filter!(x -> x ∈ elli, points)
    return points
end

include("ellipsoid_inclusion.jl")
include("ellipsoid_intersection.jl")

# compress E1 if E1∩E2≠∅
# return nothing if impossible
function compress_if_intersection(E1::Ellipsoid, E2::Ellipsoid)
    if is_intersecting(E1, E2)
        return scale_for_noninclusion_contact_point(E1, E2)
    else
        return E1
    end
end
