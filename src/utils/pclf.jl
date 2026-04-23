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

# Any template:
function approximate_sublevel_set(
    piece::AbstractPiece,
    γ::Float64;
    xmin = -2.0,
    xmax = 2.0,
    ymin = -2.0,
    ymax = 2.0,
    N = 300,
)
    xs = range(xmin, xmax; length = N)
    ys = range(ymin, ymax; length = N)
    vals = [piece_value(piece, [x, y]) for y in ys, x in xs]
    mask = vals .<= γ
    return xs, ys, vals, mask
end

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

# Polyhedral Lyapunov functions: Gx <= w
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
        push!(cons, LazySets.HalfSpace(gi, gamma * w[i]))
        push!(cons, LazySets.HalfSpace(-gi, gamma * w[i]))
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

"""
Compute a path-complete Lyapunov function (PCLF) with **quadratic (ellipsoidal) pieces**
for a switched linear system.

Each node `s` of the graph is associated with a quadratic Lyapunov function:

    V_s(x) = xᵀ P_s x,

where `P_s` is a symmetric positive definite matrix. The corresponding sublevel sets
are ellipsoids.

# Method
The method formulates a semidefinite feasibility problem (SDP) and performs a
bisection on γ. For each edge (u → v, σ), it enforces the Lyapunov inequality:

    A_σᵀ P_v A_σ ≤ γ² P_u,

implemented via linear matrix inequalities (LMIs):

    γ² P_u - A_σᵀ P_v A_σ - I ≽ 0.

Additional constraints ensure positive definiteness and boundedness of the matrices `P_s`.

# Arguments
- `f`: hybrid system containing the system matrices `A_σ`
- `G`: labeled directed graph defining the PCLF structure
- `optimizer`: JuMP-compatible SDP solver

# Keyword arguments
- `tol`: tolerance for bisection on γ
- `maxiter`: maximum number of iterations
- `MLF`: if true, extracts the Lyapunov matrices `P_s`

# Returns
- `PCLF`: structure containing the graph, Lyapunov pieces (ellipsoids), and JSR approximation

# Notes
- This method searches for a quadratic (ellipsoidal) Lyapunov function on each node.
- It relies on semidefinite programming (SDP), which is more expensive than LP-based
  polyhedral methods but often less conservative.
- The resulting Lyapunov function is smooth and globally defined on each node.
"""
function compute_quadratic_pieces_pclf(
    f::HybridSystems.HybridSystem,
    G::LabDigraph,
    optimizer;
    tol = 1e-5,
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

"""
Compute a path-complete Lyapunov function (PCLF) with **symmetric polyhedral pieces
having 2n faces** for a switched linear system.

Each node `s` of the graph is associated with a polyhedral Lyapunov function of the form:

    V_s(x) = max_i |(G_s x)_i| / w_s[i]

whose sublevel sets are polytopes:

    { x : -γ w_s ≤ G_s x ≤ γ w_s }.

# Method
The method constructs and solves a feasibility linear program (LP) using bisection on γ.
For each edge (u → v, σ), it enforces:

    |G_v A_σ G_u^{-1}| * w_u ≤ γ w_v,

where the absolute value is taken elementwise.

# Arguments
- `f`: hybrid system containing the system matrices `A_σ`
- `D`: labeled directed graph defining the PCLF structure
- `optimizer`: JuMP optimizer

# Keyword arguments
- `Gmats`: choice of matrices G_s (identity, Dict, or Vector)
- `tol`: tolerance for bisection on γ
- `maxiter`: maximum number of bisection iterations
- `MLF`: if true, extracts the Lyapunov pieces
- `verbose`: enable solver output
- `min_w`: lower bound to enforce strict positivity of w

# Returns
- `PCLF`: structure containing the graph, Lyapunov pieces, and JSR approximation

# Notes
- The resulting Lyapunov functions are structured and correspond to
  weighted ∞-norms in transformed coordinates.
- This approach is computationally efficient but may be conservative.
"""
function compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f::HybridSystems.HybridSystem,
    D::LabDigraph,
    optimizer;
    Gmats = :identity,
    tol = 1e-5,
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

"""
Compute a path-complete Lyapunov function (PCLF) with **general polyhedral pieces**
defined over a partition of the state space into cones.

Each node `s` is associated with a piecewise-linear Lyapunov function:

    V_s(x) = max_i |p_{s,i}ᵀ x|,

where the rows of a matrix `P_s` define the supporting hyperplanes of the polytope.

# Method
The method formulates a feasibility linear program (LP) based on:

1. Positivity constraints ensuring V_s(x) ≥ 0 on each cone
2. Dominance constraints ensuring correct piecewise structure
3. Decrease conditions along edges:

       V_v(A_σ x) ≤ ρ V_u(x)

These constraints are enforced on the extreme rays of the cones in `partitions`.

A bisection on ρ is used to approximate the joint spectral radius (JSR).

# Arguments
- `f`: hybrid system containing the system matrices `A_σ`
- `D`: labeled directed graph defining the PCLF structure
- `optimizer`: JuMP optimizer
- `partitions`: dictionary mapping each node to a list of cones (matrices of rays)

# Keyword arguments
- `tol`: tolerance for bisection on ρ
- `maxiter`: maximum number of iterations
- `MLF`: if true, extracts the Lyapunov pieces
- `verbose`: enable solver output
- `min_c`: lower bound on auxiliary scalar variables

# Returns
- `PCLF`: structure containing the graph, Lyapunov pieces, and JSR approximation

# Notes
- This method allows for more general polyhedral Lyapunov functions than the
  symmetric 2n-face construction.
- The number of faces depends on the number of rows of `P_s`.
- Less conservative but computationally more expensive.
- The quality depends on the chosen cone partition.
"""
function compute_polyhedral_pieces_pclf(
    f::HybridSystems.HybridSystem,
    D::LabDigraph,
    optimizer,
    partitions;
    tol = 1e-5,
    maxiter = 100,
    MLF = false,
    verbose = false,
    min_c = 1e-3,
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

    # --- vertices and indexing ---
    verts = collect(D.verts)
    l_s = length(verts)
    index_of = Dict{Any, Int}()
    for (i, v) in enumerate(verts)
        index_of[v] = i
    end

    # --- check partitions and infer dimension ---
    @assert haskey(partitions, verts[1]) "Partitions missing for vertex $(verts[1])"
    @assert !isempty(partitions[verts[1]]) "Node $(verts[1]) must have at least one cone"

    n = size(partitions[verts[1]][1], 1)

    for v in verts
        @assert haskey(partitions, v) "Partitions missing for vertex $v"
        @assert !isempty(partitions[v]) "Node $v must have at least one cone"
        for Xi in partitions[v]
            @assert size(Xi, 1) == n "All cones must live in R^n"
        end
    end

    # --- linear form helper ---
    linrow(P, i, x) = sum(P[i, k] * x[k] for k in 1:n)

    # --- solve feasibility LP for a fixed rho ---
    function feasibility_at_rho(rho::Float64; extract_solution::Bool = false)
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        # variables: one matrix P_s per node, one c_s per node
        Pvars = Dict{Any, Matrix{JuMP.VariableRef}}()
        cvars = Dict{Any, JuMP.VariableRef}()

        for v in verts
            l_v = length(partitions[v])
            Pvars[v] = JuMP.@variable(model, [1:l_v, 1:n], base_name = "P_$(index_of[v])")
            cvars[v] =
                JuMP.@variable(model, base_name = "c_$(index_of[v])", lower_bound = min_c)
        end

        # --- node-wise constraints (Theorem-style conditions) ---
        for s in verts
            Ps = Pvars[s]
            cs = cvars[s]
            cones_s = partitions[s]
            l_s_local = length(cones_s)

            # (a) positivity
            for i in 1:l_s_local
                Xi = cones_s[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    for j in 1:n
                        JuMP.@constraint(model, linrow(Ps, i, x) + cs * x[j] >= 0)
                        JuMP.@constraint(model, linrow(Ps, i, x) - cs * x[j] >= 0)
                    end
                end
            end

            # (b) dominance of row i on cone i
            for i in 1:l_s_local
                Xi = cones_s[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    for k in 1:l_s_local
                        k == i && continue
                        JuMP.@constraint(model, linrow(Ps, i, x) + linrow(Ps, k, x) >= 0)
                        JuMP.@constraint(model, linrow(Ps, i, x) - linrow(Ps, k, x) >= 0)
                    end
                end
            end
        end

        # (c) edge inequalities (s,m,d): V_d(A_m x) <= rho * V_s(x)
        for (u, v, label) in D.edges
            σ = Int(label)
            Am = A[σ]

            Pu = Pvars[u]
            Pv = Pvars[v]
            cones_u = partitions[u]
            l_u_local = length(cones_u)
            l_v_local = length(partitions[v])

            for i in 1:l_u_local
                Xi = cones_u[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    Ax = Am * x
                    for r in 1:l_v_local
                        JuMP.@constraint(
                            model,
                            rho * linrow(Pu, i, x) + linrow(Pv, r, Ax) >= 0
                        )
                        JuMP.@constraint(
                            model,
                            rho * linrow(Pu, i, x) - linrow(Pv, r, Ax) >= 0
                        )
                    end
                end
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        feasible = (st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT)

        if !feasible
            return (feasible = false, P_by_node = nothing, c_by_node = nothing)
        end

        if !extract_solution
            return (feasible = true, P_by_node = nothing, c_by_node = nothing)
        end

        P_by_node = Dict{Any, Matrix{Float64}}()
        c_by_node = Dict{Any, Float64}()

        for v in verts
            P_by_node[v] = JuMP.value.(Pvars[v])
            c_by_node[v] = JuMP.value(cvars[v])
        end

        return (feasible = true, P_by_node = P_by_node, c_by_node = c_by_node)
    end

    a = 0.0
    b = 0.0
    for Ai in A
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(Ai))))
        b = max(b, LinearAlgebra.opnorm(Ai, Inf))   # infinity norm (max row sum)
    end

    # --- bisection over rho ---
    iter = 0
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        rho_trial = (a + b) / 2
        res = feasibility_at_rho(rho_trial; extract_solution = false)

        if res.feasible
            b = rho_trial
        else
            a = rho_trial
        end
    end

    gamma = b

    # --- final solve to extract pieces ---
    pieces = Dict{Any, AbstractPiece}()
    if MLF
        final_res = feasibility_at_rho(gamma; extract_solution = true)
        if !(final_res.feasible)
            @warn "Final LP not feasible/optimal; status = infeasible"
        else
            for v in verts
                P = final_res.P_by_node[v]
                pieces[v] = PolyhedralPiece(P, ones(size(P, 1)))
            end
        end
    end

    return PCLF(D, pieces, gamma)
end

# Evaluation of a piece at a vector x:
piece_value(::AbstractPiece, ::AbstractVector{<:Real}) =
    error("piece_value not implemented")

function piece_value(p::EllipsoidalPiece, x::AbstractVector{<:Real})
    return LinearAlgebra.dot(x, p.P * x)
end

function piece_value(p::PolyhedralPiece, x::AbstractVector{<:Real})
    gx = p.G * x
    return maximum(abs.(gx) ./ p.w)
end

struct ObserverCLFPiece{U} <: AbstractPiece
    observer_states::Vector{Set{U}}
    base_pieces::Dict{U, AbstractPiece}
end

function get_sublevel_set(piece::ObserverCLFPiece, γ::Float64)
    parts = LazySets.HPolytope[]

    for S in piece.observer_states
        isempty(S) && continue

        cons = LazySets.HalfSpace[]
        for i in S
            Pi = piece.base_pieces[i]
            @assert Pi isa PolyhedralPiece

            for k in 1:size(Pi.G, 1)
                gk = vec(Pi.G[k, :])
                push!(cons, LazySets.HalfSpace(gk, γ * Pi.w[k]))
                push!(cons, LazySets.HalfSpace(-gk, γ * Pi.w[k]))
            end
        end

        push!(parts, LazySets.HPolytope(cons))
    end

    return UT.SemiLinearSet(parts)
end

function piece_value(p::ObserverCLFPiece, x::AbstractVector{<:Real})
    best = Inf
    for S in p.observer_states
        isempty(S) && continue
        worst = -Inf
        for i in S
            worst = max(worst, piece_value(p.base_pieces[i], x))
        end
        best = min(best, worst)
    end
    return best
end

graph_labels(g::LabDigraph{T, U}) where {T, U} = unique(last.(g.edges))

function successor_subset(g::LabDigraph{T, U}, S::Set{U}, h::T) where {T, U}
    Tset = Set{U}()
    for (src, dst, lab) in g.edges
        if lab == h && src in S
            push!(Tset, dst)
        end
    end
    return Tset
end

canonical_state(S::Set{U}) where {U} = Tuple(sort(collect(S); by = x -> string(x)))

function build_observer_graph(g::LabDigraph{T, U}) where {T <: Real, U}
    alphabet = collect(graph_labels(g))
    start = Set(g.verts)

    states = Vector{Set{U}}()
    trans = Dict{Tuple{Int, T}, Int}()
    seen = Dict{Any, Int}()

    push!(states, start)
    seen[canonical_state(start)] = 1

    queue = [1]
    while !isempty(queue)
        k = popfirst!(queue)
        S = states[k]

        for h in alphabet
            Tset = successor_subset(g, S, h)
            isempty(Tset) && continue

            key = canonical_state(Tset)
            if !haskey(seen, key)
                seen[key] = length(states) + 1
                push!(states, Tset)
                push!(queue, length(states))
            end
            trans[(k, h)] = seen[key]
        end
    end

    return states, trans, alphabet
end

function common_lyapunov_graph(labels::Vector{T}) where {T <: Real}
    node = :clf
    verts = Set([node])
    edges = [(node, node, h) for h in labels]
    return LabDigraph{T, Symbol}(edges, verts)
end

function build_common_lyapunov(pclf::PCLF)
    states, _, alphabet = build_observer_graph(pclf.graph)

    U = eltype(pclf.graph.verts)
    pieces_typed = Dict{U, AbstractPiece}(pclf.pieces)

    clf_piece = ObserverCLFPiece(states, pieces_typed)
    clf_graph = common_lyapunov_graph(alphabet)

    return PCLF(clf_graph, Dict(:clf => clf_piece), pclf.JSRapprox)
end

end # module
