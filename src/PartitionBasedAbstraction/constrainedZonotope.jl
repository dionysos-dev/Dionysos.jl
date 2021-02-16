import LazySets.is_intersection_empty,LazySets.intersection,LazySets.center
import RecipesBase.plot!

mutable struct ConstrainedZonotope{N} <: AbstractPolytope{N}
    c::Vector{N}
    G::Array{N,2}
    A::Array{N,2}
    b::Vector{N}
    P::HPolytope
end
ConstrainedZonotope(c,G,A,b) =  ConstrainedZonotope{Float64}(c,G,A,b,convert_CG_H(c,G,A,b))

# H-rep → CG-rep or V-rep → CG-rep
function convert_HorV_CG(P::AbstractPolytope)

    # calculate the vertices of the polytope
    v = vertices_list(P)
    n = length(v[1])
    v_m = [v[j][i] for i=1:n, j=1:length(v)]
    # calculate a bounding box for the constrained zonotope and convert it to a zonotope
    box = Hyperrectangle(low=[minimum(v_m[i,:]) for i=1:n], high=[maximum(v_m[i,:]) for i=1:n])
    Z = convert(Zonotope,box)
    c = LazySets.center(Z)
    G = genmat(Z)
    # extract H-rep data of p
    D = constraints_list(P)
    A = [D[i].a[j] for i=1:length(D), j=1:n]
    b = [D[i].b for i=1:length(D)]
    # calculate the lower bound sigma for A*x \in [sigma,b]
    M = A*v_m
    σ = [minimum(M[i,:]) for i=1:length(D)]
    # construct constrained zonotope object according to equation
    G_ = hcat(G, zeros(size(G,1),size(A,1)))
    A_ = hcat(A*G, diagm((σ-b)./2))
    b_ = (b+σ)./2 - A*c

    return ConstrainedZonotope(c,G_,A_,b_)
end

# CG-rep → H-rep
function convert_CG_H(CZ::ConstrainedZonotope)
    return convert_CG_H(CZ.c,CZ.G,CZ.A,CZ.b)
end

function convert_CG_H(c::Vector{N},G::Array{N,2},A::Array{N,2},b::Vector{N}) where N
    n = length(c)
    g = size(G,2)

    #compute the lifted zonotope
    c_l = vcat(c,-b)
    G_l = vcat(G,A)
    Z_l = Zonotope(c_l,G_l)
    #compute the H-rep of the lifted zonotope
    HPol_l = convert(HPolytope,Z_l)
    D = constraints_list(HPol_l)
    #compute the H-rep of the constrained zonotope
    H = [D[i].a[j] for i=1:length(D), j=1:n]
    e = [D[i].b for i=1:length(D)]
    P = HPolytope(H,e)
    #eliminate the redondancy constraints
    remove_redundant_constraints!(P)

    return P
end

function LazySets.is_intersection_empty(P::AbstractPolytope,CZ::ConstrainedZonotope)
    return is_intersection_empty(P,CZ.P)
end

function LazySets.is_intersection_empty(CZ1::ConstrainedZonotope,CZ2::ConstrainedZonotope)
    return is_intersection_empty(CZ1.P,CZ2.P)
end

function LazySets.intersection(P::AbstractPolytope,CZ::ConstrainedZonotope)
    return intersection(P,CZ.P)
end

function LazySets.intersection(CZ1::ConstrainedZonotope,CZ2::ConstrainedZonotope)
    return intersection(CZ1.P,CZ2.P)
end

function LazySets.issubset(CZ::ConstrainedZonotope,P::AbstractPolytope)
    return issubset(CZ.P,P)
end

function LazySets.issubset(P::AbstractHyperrectangle,CZ::ConstrainedZonotope)
    return issubset(P,CZ.P)
end
function LazySets.issubset(P::Singleton,CZ::ConstrainedZonotope)
    return issubset(P,CZ.P)
end

function LazySets.issubset(CZ1::ConstrainedZonotope,CZ2::ConstrainedZonotope)
    return issubset(CZ1.P,CZ2.P)
end

function LazySets.center(CZ::ConstrainedZonotope)
    return CZ.c
end

function LazySets.in(x::AbstractVector{N},CZ::ConstrainedZonotope{N}) where N
    return x∈CZ.P
end

function LazySets.overapproximate(CZ::ConstrainedZonotope)
    return overapproximate(CZ.P)
end

function plot!(CZ::ConstrainedZonotope)
    plot!(CZ.P)
end

function expansion(α::Float64,Z::Zonotope)
    return Zonotope(Z.center,α.*Z.generators)
end
function expansion(α::Float64,Z::ConstrainedZonotope)
    return ConstrainedZonotope(Z.c,α.*Z.G,Z.A,Z.b)
end
function expansion(α::Float64,P::HPolytope)
    CZ = convert_HorV_CG(P)
    return convert_CG_H(expansion(α,CZ))
end

function build_generators(i,centers,neigh)
    n = length(centers[1])
    g = length(neigh)
    @assert g >= n "Number of generator must be greater or" *
            " equal than the dimension of the state space"
    G = zeros((n,g))
    for j=1:length(neigh)
        G[:,j] = 0.5*(centers[neigh[j]]-centers[i])
    end
    return G
end
