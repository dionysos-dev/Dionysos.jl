module PathCompleteFramework

import Dionysos
const UT = Dionysos.Utils
import HybridSystems
import JuMP
import MathOptInterface
const MOI = MathOptInterface
import LinearAlgebra

import LazySets

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
    G::Matrix{Float64}      # m x n matrix of rank n
    w::Vector{Float64}      # n-dimensional positive vector
end

# P(γ) = { x : -γ w <= G x <= γ w }.
function get_sublevel_set(piece::PolyhedralPiece, gamma::Float64)
    G = piece.G
    w = piece.w

    m, n = size(G)
    @assert length(w) == m "w must have length equal to number of rows of G"

    cons = LazySets.HalfSpace[]
    for i in 1:m
        gi = vec(G[i, :])
        push!(cons, LazySets.HalfSpace( gi,  gamma * w[i]))
        push!(cons, LazySets.HalfSpace(-gi,  gamma * w[i]))
    end
    return LazySets.HPolytope(cons)
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

function compute_polyhedral_pieces_pclf(
    f::HybridSystems.HybridSystem,
    D::LabDigraph,
    optimizer;
    Gmats = :identity,
    tol = 1e-8,
    maxiter = 100,
    MLF = false,
    verbose = false,
    min_w = 1e-3,
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

    # --- vertices and indexing for w variables ---
    verts = collect(D.verts)
    l_s = length(verts)
    index_of = Dict{typeof(verts[1]), Int}()
    for (i, v) in enumerate(verts)
        index_of[v] = i
    end

    # Normalize Gmats input to a Vector indexed 1..l_s
    # --- Build G_by_idx ---
    n = size(A[1], 1)
    G_by_idx = Vector{Matrix{Float64}}(undef, l_s)

    if Gmats === :identity
        # Default: all G_s = I
        for i in 1:l_s
            G_by_idx[i] = Matrix{Float64}(LinearAlgebra.I, n, n)
        end

    elseif isa(Gmats, Dict)
        for (i, v) in enumerate(verts)
            @assert haskey(Gmats, v) "Gmats missing vertex $v"
            G_by_idx[i] = Array(Gmats[v])
        end

    elseif isa(Gmats, AbstractVector)
        @assert length(Gmats) >= l_s "Gmats vector too short"
        for i in 1:l_s
            G_by_idx[i] = Array(Gmats[i])
        end

    else
        error("Gmats must be :identity, Dict or Vector of matrices")
    end

    # sizes and inverses
    n = size(G_by_idx[1], 1)
    for i in 1:l_s
        @assert size(G_by_idx[i], 1) == n && size(G_by_idx[i], 2) == n "All G matrices must be n×n"
    end
    Ginv = [inv(G_by_idx[i]) for i in 1:l_s]

    # --- Precompute M matrices for each edge: M = |G_d * A_sigma * G_s^{-1}| ---
    edge_list = D.edges  # expected iterable of (u,v,label)
    Mlist = Vector{Tuple{Int, Int, Int, Matrix{Float64}}}()  # (ui,vi,sigma,M)
    for (u, v, label) in edge_list
        ui = index_of[u];
        vi = index_of[v];
        σ = Int(label)
        M = abs.(G_by_idx[vi] * A[σ] * Ginv[ui])   # nonnegative n×n matrix
        push!(Mlist, (ui, vi, σ, M))
    end

    # --- initial upper bound b: max row-sum among M matrices (finite) ---
    a = 0.0
    b = 0.0
    for Ai in A
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(Ai))))
        b = max(b, LinearAlgebra.opnorm(Ai, Inf))   # infinity norm (max row sum)
    end

    # --- bisection ---
    iter = 0
    feasible_at = false
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        gamma = (a + b) / 2

        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        # variables: w_i in R^n with strict positivity lower bound min_w
        wvars = [JuMP.@variable(model, [1:n]) for i in 1:l_s]

        # enforce strict positivity
        for i in 1:l_s, k in 1:n
            JuMP.@constraint(model, wvars[i][k] >= min_w)
        end

        # constraints: for every precomputed M (ui->vi): M * w_ui <= gamma * w_vi
        for (ui, vi, σ, M) in Mlist
            for k in 1:n
                JuMP.@constraint(
                    model,
                    sum(M[k, p] * wvars[ui][p] for p in 1:n) <= gamma * wvars[vi][k]
                )
            end
        end

        # Feasibility check (no objective)
        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        feasible_at = (st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT)

        if feasible_at
            b = gamma
        else
            a = gamma
        end
    end

    gamma = b

    # --- final solve to extract pieces (if requested) ---
    pieces = Dict{Any, AbstractPiece}()
    if MLF
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        wvars = [JuMP.@variable(model, [1:n]) for i in 1:l_s]

        # enforce strict positivity
        for i in 1:l_s, k in 1:n
            JuMP.@constraint(model, wvars[i][k] >= min_w)
        end

        for (ui, vi, σ, M) in Mlist
            for k in 1:n
                JuMP.@constraint(
                    model,
                    sum(M[k, p] * wvars[ui][p] for p in 1:n) <= gamma * wvars[vi][k]
                )
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        if !(st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT)
            @warn "Final polyhedral LP not feasible/optimal; status = $st"
        else
            for (i, v) in enumerate(verts)
                wnum = JuMP.value.(wvars[i])
                # numeric safety: enforce positivity
                wnum .= max.(wnum, min_w)
                pieces[v] = PolyhedralPiece(G_by_idx[i], Array(wnum))
            end
        end
    end

    return PCLF(D, pieces, gamma)
end

end # module
