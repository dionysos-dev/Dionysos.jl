module PathCompleteFramework

import Dionysos
const UT = Dionysos.Utils
import HybridSystems
import JuMP
import MathOptInterface
const MOI = MathOptInterface
import LinearAlgebra

"""
    struct LabDigraph{T<:Real, U}

Store a graph as an explicit list of edges (u, v, label),
preserving parallel edges and arbitrary vertex types.
"""
struct LabDigraph{T <: Real, U}
    edges::Vector{Tuple{U, U, T}}
    verts::Set{U}
end

function edgeList_to_LabDigraph(edges::Vector{Tuple{U, U, T}}) where {T <: Real, U}
    verts = Set{U}(v for e in edges for v in (e[1], e[2]))
    return LabDigraph{T, U}(edges, verts)
end

abstract type AbstractPiece end

function get_sublevel_set(piece::AbstractPiece, gamma::Float64) end

# Quadratic Lyapunov functions:
mutable struct EllipsoidalPiece <: AbstractPiece
    P::Matrix{Float64}   # symmetric positive-definite matrix
end

function get_sublevel_set(piece::EllipsoidalPiece, gamma::Float64)
    n = size(piece.P, 1)
    center = zeros(n)
    elli = UT.Ellipsoid(piece.P, center)
    return UT.get_sublevel_set(elli, gamma)
end

# Polyhedral Lyapunov functions:
mutable struct PolyhedralPiece <: AbstractPiece
    h::Any
end

function get_sublevel_set(piece::PolyhedralPiece, gamma::Float64)
    return
end

"""
    mutable struct PCLF

Store a path-complete Lyapunov function (i.e. a graph and a collection of 
Lyapunov pieces) for a linear switched system and the corresponding 
JSR approximation
"""
mutable struct PCLF
    graph::LabDigraph
    pieces::Dict{Any, AbstractPiece}
    JSRapprox::Float64
end

function generate_DeBruijn_edges(M::Int, k::Int; dual::Bool = false)
    @assert M ≥ 1
    @assert k ≥ 0

    # Common Lyapunov function graph:
    if k == 0
        edges = Vector{Tuple{Int, Int, Int}}()
        v = (1)
        for s in 1:M
            push!(edges, (v, v, s))
        end
        return edgeList_to_LabDigraph(edges)
    end

    # Otherwise: 
    edges = Vector{Tuple{NTuple{k, Int}, NTuple{k, Int}, Int}}()

    iterables = ntuple(_ -> 1:M, k)
    nodes = collect(Iterators.product(iterables...))

    for node in nodes
        for state in nodes
            if node[2:end] == state[1:(end - 1)]
                if !dual
                    push!(edges, (node, state, state[end]))
                else
                    push!(edges, (state, node, state[end]))
                end
            end
        end
    end

    return edgeList_to_LabDigraph(edges)
end

function compute_quadratic_pieces_pclf(
    f::HybridSystems.HybridSystem,
    G::LabDigraph,
    optimizer;
    tol = 1e-8,
    maxiter = 200,
    MLF = false,
)

    # --- extract matrices from resetmaps ---
    RMs = f.resetmaps
    A = Vector{Matrix{Float64}}(undef, length(RMs))
    for (i, rm) in enumerate(RMs)
        if isa(rm, AbstractMatrix)
            A[i] = Array(rm)
        elseif :A in fieldnames(typeof(rm))
            A[i] = Array(getfield(rm, :A))
        else
            error("Cannot extract matrix from resetmap of type $(typeof(rm)).")
        end
    end

    # --- vertices and indexing for P variables ---
    verts = collect(G.verts)
    l_s = length(verts)
    index_of = Dict{typeof(verts[1]), Int}()
    for (i, v) in enumerate(verts)
        index_of[v] = i
    end

    # --- edge list from LabDigraph ---
    edge_list = G.edges

    # --- matrix dimension checks / Identity ---
    n = size(A[1], 1)
    I_n = Matrix{Float64}(LinearAlgebra.I, n, n)

    # --- initial bounds a and b ---
    a = 0.0
    b = 0.0
    for Ai in A
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(Ai))))
        b = max(b, LinearAlgebra.opnorm(Ai, 2))
    end

    # --- bisection loop ---
    iter = 0
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        gamma = (a + b) / 2
        γ2 = gamma^2

        model = JuMP.Model(optimizer)

        # create P variables (anonymous) as matrix variables
        P = [JuMP.@variable(model, [1:n, 1:n], Symmetric) for i in 1:l_s]

        # add initial lower bound: P[i] >= 0.5*I
        for i in 1:l_s
            JuMP.@constraint(model, P[i] - 0.5 * I_n in JuMP.PSDCone())
            JuMP.@constraint(model, 1000*I_n - P[i] in JuMP.PSDCone())
        end

        # add LMIs for each edge (u -> v with label)
        for (u, v, label) in edge_list
            Pu_idx = index_of[u]
            Pv_idx = index_of[v]

            # assume label is integer index into A; adapt if labels are different
            Albl = A[Int(label)]
            expr = γ2 * P[Pu_idx] - (Albl' * P[Pv_idx] * Albl) - I_n
            JuMP.@constraint(model, expr in JuMP.PSDCone())
        end

        JuMP.set_silent(model)
        JuMP.optimize!(model)

        st = JuMP.termination_status(model)
        if st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT
            b = gamma
        else
            a = gamma
        end
    end

    gamma = b

    # --- final solve to extract P (if requested) ---
    pieces = Dict{typeof(verts[1]), AbstractPiece}()

    if MLF
        model = JuMP.Model(optimizer)

        # create vars again for final solve
        P = [JuMP.@variable(model, [1:n, 1:n], Symmetric) for i in 1:l_s]
        for i in 1:l_s
            JuMP.@constraint(model, P[i] - 0.5 * I_n in JuMP.PSDCone())
        end

        γ2 = gamma^2
        for (u, v, label) in edge_list
            Pu_idx = index_of[u]
            Pv_idx = index_of[v]
            Albl = A[Int(label)]
            expr = γ2 * P[Pu_idx] - (Albl' * P[Pv_idx] * Albl) - I_n
            JuMP.@constraint(model, expr in JuMP.PSDCone())
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        if st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT
            # build EllipsoidalPiece objects keyed by original vertex ids
            for (i, v) in enumerate(verts)
                Pnum = JuMP.value.(P[i])
                # Optionally enforce symmetry numerically: symmetrize small numerical noise
                Pnum = 0.5 * (Pnum + Pnum')
                pieces[v] = EllipsoidalPiece(Pnum)
            end
        else
            @warn "Final solve not feasible/optimal. status = $st"
            pieces = Dict{typeof(verts[1]), AbstractPiece}()
        end
    end

    return PCLF(G, pieces, gamma)
end

function compute_polyhedral_pieces_pclf(f::HybridSystems.HybridSystem, G::LabDigraph)
    return
end

end # module
