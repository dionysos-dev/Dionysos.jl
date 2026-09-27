module PathCompleteFramework

import ..Utils
const UT = Utils
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

"""
    is_complete(G::LabDigraph, modes) -> Bool

Whether every node of `G` has an outgoing edge for every mode in `modes`.

Completeness and path-completeness are different properties and license different uses. A
path-complete graph represents every finite mode *word* by some path, which is what synthesis needs:
a co-safe specification is witnessed in finite time, and the witness word has a path to follow.
Verification under arbitrary switching needs more — a node missing an outgoing edge for some mode
removes a move the real adversary has, so the abstraction can certify a universal property the
concrete system does not satisfy. Completeness is exactly what closes that gap.

See also [`is_deterministic`](@ref), which bounds the same counts from above rather than below.
"""
function is_complete(G::LabDigraph, modes)
    return all(any(e -> e[1] == s && e[3] == m, G.edges) for s in G.verts for m in modes)
end

"""
    is_deterministic(G::LabDigraph, modes) -> Bool

Whether every `(node, mode)` pair of `G` has at most one outgoing edge.

When a graph branches, the successor node is not determined by the mode, and something must resolve
the choice. A controller that tracks one node resolves it in its own favour; an adversary that
resolves it instead removes freedom the concrete system grants. The distinction matters for
synthesis guarantees, so it is worth being able to test.

See also [`is_complete`](@ref).
"""
function is_deterministic(G::LabDigraph, modes)
    return all(
        count(e -> e[1] == s && e[3] == m, G.edges) <= 1 for s in G.verts for m in modes
    )
end

"""
    is_co_complete(G::LabDigraph, modes) -> Bool

Whether every node of `G` has an **incoming** edge for every mode in `modes`.

The mirror of [`is_complete`](@ref), and the reason both are worth testing: each licenses a different
closed-form common Lyapunov function, and the two are genuinely different objects.

Write the edge condition as `V_d(A_m x) ≤ γ V_s(x)` for every edge `(s, m, d)`.

- **complete** ⟹ `V_min = min_i V_i` is a common Lyapunov function. If the minimum at `x` is attained
  at `i`, completeness supplies an edge `(i, m, d)`, so `min_j V_j(A_m x) ≤ V_d(A_m x) ≤ γ V_i(x)`.
  Its sublevel set is the **union** of the pieces' sublevel sets — non-convex, but the larger set.
- **co-complete** ⟹ `V_max = max_i V_i` is one. For each target `i`, co-completeness supplies an edge
  `(s, m, i)`, so `V_i(A_m x) ≤ γ V_s(x) ≤ γ V_max(x)`. Its sublevel set is the **intersection** —
  convex and cheap, but the smaller set.

On De Bruijn graphs exactly one holds: the primal (a node records the last mode played) is complete
and not co-complete; the dual (a node commits to the mode played next, and any node may follow) is
co-complete and not complete. So determinising a dual graph necessarily *loses coverage*, and
determinising a primal one necessarily *loses convexity*. See [`build_common_lyapunov`](@ref).
"""
function is_co_complete(G::LabDigraph, modes)
    return all(any(e -> e[2] == s && e[3] == m, G.edges) for s in G.verts for m in modes)
end

"""
    enabled_modes(G::LabDigraph) -> Dict{U, Set{Int}}

The modes each node of `G` may actually emit.

This is the language talking, and it has to be carried forward explicitly rather than rediscovered
later from the absence of a transition. Downstream, "node `q` has no successor under mode `m`" has two
irreconcilable causes: the language forbids `m` at `q`, so the environment *cannot* play it; or `m` is
allowed but its successor left the region the abstraction covers, so a real behaviour was dropped. A
universal (`∀`) answer must ignore the first and pessimistically account for the second, and once the
graph is gone the two are indistinguishable. See [`is_complete`](@ref) for when the distinction is
vacuous: on a complete graph every node enables every mode and nothing is ever forbidden.
"""
function enabled_modes(G::LabDigraph{T, U}) where {T, U}
    out = Dict{U, Set{Int}}(v => Set{Int}() for v in G.verts)
    for (u, _, label) in G.edges
        push!(out[u], Int(label))
    end
    return out
end

"""
    restricts_future(G::LabDigraph) -> Bool

Whether the node reached so far changes which modes are available next.

`false` on a De Bruijn graph, where a node records history only and every node enables every mode, so
a result holding at one node holds at all of them. `true` for a genuine language restriction such as
"no two consecutive uses of mode 1", where the node constrains the *future* — and then results at
different nodes are answers to different questions and must not be merged without deciding which
node a run may start in.
"""
function restricts_future(G::LabDigraph)
    enabled = enabled_modes(G)
    isempty(enabled) && return false
    reference = first(values(enabled))
    return any(!=(reference), values(enabled))
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
    # The PCLF piece is P-natural; one inversion at construction.
    elli = LazySets.Ellipsoid(center, UT._symmetrize(inv(piece.P)))
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
    A = UT.mode_matrices(f)

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
        b = max(b, LinearAlgebra.opnorm(Ai, 2))
    end

    # The lower bracket must bound the rate over the *graph's* language, not the system's. A mode
    # carrying a self-loop can be repeated for ever, so the language contains `mᵂ` and the rate is
    # at least `ρ(A_m)`; a mode with no self-loop may be forbidden from repeating, and assuming
    # otherwise clamps the bisection above the true answer -- silently, since it then converges to
    # the bracket floor for every template. With no self-loop this degrades to 0, and for a single
    # node carrying every mode it reproduces the previous bracket exactly.
    for (u, v, label) in G.edges
        u == v || continue
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(A[Int(label)]))))
    end

    # --- bisection loop ---
    iter = 0
    any_feasible = false
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
            any_feasible = true
            b = gamma
        else
            a = gamma
        end
    end

    # See the note in `compute_polyhedral_pieces_pclf`: `b` is only an unverified norm bound until
    # some trial certifies it, so returning it unconditionally reports a failure as a result.
    if !any_feasible
        @warn "compute_quadratic_pieces_pclf: no feasible contraction rate found; no certificate \
               exists for this graph and template. Returning JSRapprox = Inf."
        return PCLF(G, Dict{typeof(verts[1]), AbstractPiece}(), Inf)
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

# Reference
This symmetric `2n`-face construction of polyhedral path-complete Lyapunov
functions follows [athanasopoulos2019polyhedral](@cite).
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
    A = UT.mode_matrices(f)

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
        ui = index_of[u]
        vi = index_of[v]
        σ = Int(label)
        M = abs.(G_by_idx[vi] * A[σ] * Ginv[ui])   # nonnegative n×n matrix
        push!(Mlist, (ui, vi, σ, M))
    end

    # --- initial upper bound b: max row sum among the M matrices ---
    #
    # It must be taken on the `M = |G_v A_σ G_u^{-1}|`, NOT on the raw `A_σ`. The constraint the
    # bisection tests is `M w_u <= γ w_v`, so giving every `w` the same value makes it feasible
    # exactly at `γ = max row sum of M` — a bracket the search is *guaranteed* to start inside.
    # `opnorm(A_σ, Inf)` is that same quantity only when the templates are the identity; under a
    # rotated template it can be strictly smaller, and then the top of the bracket is itself
    # infeasible and a certifiable system is reported as having no certificate.
    a = 0.0
    b = 0.0
    for (_, _, _, M) in Mlist
        b = max(b, maximum(sum(M; dims = 2)))
    end

    # See `compute_quadratic_pieces_pclf` for why the lower bracket is language-aware.
    for (u, v, label) in D.edges
        u == v || continue
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(A[Int(label)]))))
    end

    # Whether one trial rate admits piece weights, as a closure so the bracket can be tested at `b`
    # itself rather than only at points strictly inside it.
    function feasible_at_rate(gamma)
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        # variables: w_i in R^n with strict positivity lower bound min_w
        wvars = [JuMP.@variable(model, [1:n]) for i in 1:l_s]
        for i in 1:l_s, k in 1:n
            JuMP.@constraint(model, wvars[i][k] >= min_w)
        end

        # constraints: for every precomputed M (ui->vi): M * w_ui <= gamma * w_vi
        for (ui, vi, _, M) in Mlist
            for k in 1:n
                JuMP.@constraint(
                    model,
                    sum(M[k, p] * wvars[ui][p] for p in 1:n) <= gamma * wvars[vi][k]
                )
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        return st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT
    end

    # --- bisection ---
    #
    # The top of the bracket is TESTED, not assumed feasible. Without this the search reports failure
    # whenever the optimum sits exactly at `b`, because the lower bracket `a = ρ(A_m)` over self-loops
    # then equals `b`, `b - a > tol` is false at once, the loop body never runs and `any_feasible`
    # stays false. That is not a corner case: `A = diag(0.6, 0.48)` under a π/6-rotated template has
    # Perron root exactly 0.6 — a similarity transform preserves the spectrum — so `a` and `b` meet
    # on the true answer and a perfectly good certificate was reported as absent.
    any_feasible = feasible_at_rate(b)

    iter = 0
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        gamma = (a + b) / 2
        if feasible_at_rate(gamma)
            any_feasible = true
            b = gamma
        else
            a = gamma
        end
    end

    # With a fixed template and free positive weights a feasible rate ALWAYS exists — give every `w`
    # the same value and the constraints hold at `γ = max row sum of M`, which is how `b` is chosen.
    # So reaching here with nothing feasible does not mean "no certificate"; it means the solver could
    # not certify even the bracket that is feasible by construction, which is a numerical failure and
    # is reported as such. `Inf` is kept for that case rather than returning a rate no solve backed.
    if !any_feasible
        @warn "compute_symmetric_2n_faces_polyhedral_pieces_pclf: the solver rejected even the \
               guaranteed bracket γ = $(b), which is feasible by construction (equal weights). This \
               is a numerical failure, not an absent certificate. Returning JSRapprox = Inf."
        return PCLF(D, Dict{Any, AbstractPiece}(), Inf)
    end

    gamma = b

    # A rate is always returned, and it is the caller's job to ask whether it contracts. Saying so
    # here because the previous version conflated "no contraction" with "no answer", and a rate of,
    # say, 1.2 is strictly more useful than `Inf`: it says how far the template is from working.
    if gamma >= 1
        @warn "compute_symmetric_2n_faces_polyhedral_pieces_pclf: best rate for this graph and \
               template is $(gamma) >= 1, so it certifies no contraction." maxlog = 1
    end

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
    A = UT.mode_matrices(f)

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
        b = max(b, LinearAlgebra.opnorm(Ai, Inf))   # infinity norm (max row sum)
    end

    # See `compute_quadratic_pieces_pclf` for why the lower bracket is language-aware.
    for (u, v, label) in D.edges
        u == v || continue
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(A[Int(label)]))))
    end

    # --- bisection over rho ---
    # `b` starts at a trivial norm bound that is *not* known to be feasible. If no trial rate is
    # ever feasible, `b` never moves, and returning it would report that trivial bound as though it
    # were a computed certificate — silently, and with the same value for every template, which is
    # how this failure previously masqueraded as a non-monotone bound across template refinements.
    iter = 0
    any_feasible = false
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        rho_trial = (a + b) / 2
        res = feasibility_at_rho(rho_trial; extract_solution = false)

        if res.feasible
            any_feasible = true
            b = rho_trial
        else
            a = rho_trial
        end
    end

    if !any_feasible
        @warn "compute_polyhedral_pieces_pclf: no feasible contraction rate found; no certificate \
               exists for this graph and template. Returning JSRapprox = Inf."
        return PCLF(D, Dict{Any, AbstractPiece}(), Inf)
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

function get_sublevel_set(piece::ObserverCLFPiece, γ::Float64; atol::Float64 = 1e-6)
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

        P = UT.clean_poly(LazySets.HPolytope(cons))

        if isempty(parts)
            push!(parts, P)
        else
            prev_union = UT.semilinear_set(parts)
            remainder = UT.set_difference_decompose(P, prev_union; atol = atol)
            append!(parts, remainder)
        end
    end

    return UT.semilinear_set(parts)
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

# Drop observer states that CONTAIN another one.
#
# The sublevel set is a union over states of an intersection over the nodes of each, so `S′ ⊆ S`
# makes `S` redundant: `∩_{i∈S} ⊆ ∩_{i∈S′}`. Dropping it leaves the function itself unchanged, since
# in `min_S max_{i∈S} V_i` a superset's `max` is never the smallest, but removes one polytope from
# the union, and with it the disjoint-decomposition work that polytope would have caused downstream.
#
# This is not a micro-optimisation. On a complete graph the observer reaches the full vertex set (the
# initial uncertainty) AND every singleton; the full set contributes the intersection, which then has
# to be carved out of each singleton's piece, turning two convex polytopes into five disjoint parts.
# Every cell of a quotient built on such a certificate inherits that fragmentation.
#
# A comment rather than a docstring: the `@autodocs` blocks filter underscore-prefixed internals out
# of the manual, while `checkdocs = :all` demands that every docstring appear in it, so a private
# helper cannot carry one.
function _drop_redundant_supersets(states::Vector{Set{U}}) where {U}
    keep = Set{U}[]
    for S in states
        isempty(S) && continue
        any(other -> other != S && issubset(other, S), states) && continue
        push!(keep, S)
    end
    return isempty(keep) ? states : keep
end

"""
    build_common_lyapunov(pclf::PCLF; mode = :auto) -> PCLF

The common Lyapunov function induced by `pclf`, as a one-node `PCLF`.

`mode` selects the construction, and the default checks the graph's structure rather than always
paying for the subset construction:

| `mode` | needs | `V*` | sublevel set |
| :--- | :--- | :--- | :--- |
| `:min` | [`is_complete`](@ref) | `min_i V_i` | the **union** of the pieces — non-convex, the larger set |
| `:max` | [`is_co_complete`](@ref) | `max_i V_i` | the **intersection** — convex, the smaller set |
| `:observer` | nothing | `min_S max_{i∈S} V_i` | union over the observer's reachable subsets |
| `:auto` | — | the first of the three that applies | |

`:auto` prefers `:min` when both structural tests pass, because a certificate is judged first by what
it certifies and the union is the larger region; `:max` is then available explicitly and is cheaper.
On a De Bruijn graph the question does not arise — the primal is complete and not co-complete, the
dual the reverse — so each admits exactly one closed form.

The observer fallback is not wrong on the structured graphs, merely wasteful: it recovers the same
function (on a complete graph its singletons make `min_S max_{i∈S} V_i` collapse to `min_i V_i`) but
enumerates redundant subsets on the way. Those are filtered out in every mode.
"""
function build_common_lyapunov(pclf::PCLF; mode::Symbol = :auto)
    g = pclf.graph
    labels = graph_labels(g)
    alphabet = collect(labels)

    if mode === :auto
        mode = if is_complete(g, labels)
            :min
        elseif is_co_complete(g, labels)
            :max
        else
            :observer
        end
    end

    U = eltype(g.verts)
    states = if mode === :min
        is_complete(g, labels) || error(
            "`mode = :min` needs a COMPLETE graph (every node an outgoing edge for every mode); " *
            "min_i V_i is not a common Lyapunov function otherwise.",
        )
        [Set{U}([v]) for v in g.verts]
    elseif mode === :max
        is_co_complete(g, labels) || error(
            "`mode = :max` needs a CO-COMPLETE graph (every node an incoming edge for every mode); " *
            "max_i V_i is not a common Lyapunov function otherwise.",
        )
        [Set{U}(g.verts)]
    elseif mode === :observer
        first(build_observer_graph(g))
    else
        error("unknown mode $(repr(mode)); expected :auto, :min, :max or :observer")
    end

    pieces_typed = Dict{U, AbstractPiece}(pclf.pieces)
    clf_piece = ObserverCLFPiece(_drop_redundant_supersets(states), pieces_typed)
    clf_graph = common_lyapunov_graph(alphabet)

    return PCLF(clf_graph, Dict(:clf => clf_piece), pclf.JSRapprox)
end

function base_conic_partition_2d()
    v1 = [1.0, 0.0]
    v2 = [1.0, 1.0]
    v3 = [0.0, 1.0]
    v4 = [-1.0, 1.0]
    v5 = [-1.0, 0.0]

    return [hcat(v1, v2), hcat(v2, v3), hcat(v3, v4), hcat(v4, v5)]
end

function split_cone(C::AbstractMatrix)
    a = vec(C[:, 1])
    b = vec(C[:, 2])
    m = a + b
    return hcat(a, m), hcat(m, b)
end

function conic_partitions_2d(order::Int)
    order >= 1 || error("order must be >= 1")

    cones = base_conic_partition_2d()

    for _ in 2:order
        refined = Matrix{Float64}[]
        for C in cones
            C1, C2 = split_cone(C)
            push!(refined, C1)
            push!(refined, C2)
        end
        cones = refined
    end

    return cones
end

function conic_partitions_dict_2d(order::Int, node_ids)
    cones = conic_partitions_2d(order)
    return Dict(id => cones for id in node_ids)
end

# ------------------------------------------------------------
# Checking a certificate
# ------------------------------------------------------------

# Piece values are homogeneous, but not all of the same degree: `max|Gx|/w` scales linearly in `x`
# while `x'Px` scales quadratically, and the solvers constrain them accordingly — `γ` in the
# polyhedral programs, `γ²` in the semidefinite one. `JSRapprox` is the rate of the gauge in both
# cases, so an observed value ratio must be taken to the power `1/degree` before it can be compared
# with it; without that a valid quadratic certificate reads as `ρ²` and looks better than it is.
piece_degree(::AbstractPiece) = error("piece_degree not implemented")
piece_degree(::PolyhedralPiece) = 1
piece_degree(::EllipsoidalPiece) = 2
piece_degree(p::ObserverCLFPiece) = piece_degree(first(values(p.base_pieces)))

# Both checkers report the same shape, and both compare against `JSRapprox` with a relative
# tolerance, since a rate returned by a solver is attained only up to its own accuracy.
function _rate_report(by_edge, rho, rtol)
    rate = maximum(values(by_edge))
    violated = sort!([e for (e, r) in by_edge if r > rho * (1 + rtol)]; by = string)
    return (; rate = rate, by_edge = by_edge, violated = violated, certified_rate = rho)
end

function _edge_pieces(pclf::PCLF, u, v)
    haskey(pclf.pieces, u) && haskey(pclf.pieces, v) ||
        error("The certificate has no piece for node $u or node $v.")
    pu, pv = pclf.pieces[u], pclf.pieces[v]
    piece_degree(pu) == piece_degree(pv) || error(
        "Nodes $u and $v carry pieces of different homogeneity degree; their values are not \
         comparable along an edge.",
    )
    return pu, pv
end

"""
    check_pclf(pclf, f; nsample = 20_000, rng = nothing, rtol = 1e-6)
    check_pclf(pclf, A; nsample = 20_000, rng = nothing, rtol = 1e-6)

Test the contraction inequality that defines the certificate `pclf`, by sampling.

Along every edge `(s, s′, m)` of the graph a path-complete Lyapunov function claims

    V_{s′}(A_m x) ≤ ρ^d V_s(x)    for all x,

with `ρ = pclf.JSRapprox` and `d` the homogeneity degree of the pieces. Both sides are homogeneous,
so the claim does not depend on the scale of `x` and sampling directions covers it.

The mode matrices are taken either from a `HybridSystems.HybridSystem` or given directly as a
vector indexed by mode label.

Returns `(; rate, by_edge, violated, certified_rate)`. `rate` is the largest rate observed over all
edges and `by_edge` the rate of each, both expressed on the gauge — that is, as
`(V_{s′}(A_m x) / V_s(x))^(1/d)` — so that they are directly comparable with `pclf.JSRapprox`
whatever the piece type. `violated` lists the edges exceeding it by more than `rtol` in relative
terms.

A sample underestimates a supremum, so this refutes a certificate rather than establishing one: a
`rate` above `ρ` proves the certificate wrong, a `rate` below it is evidence and not a proof. When
the pieces are polyhedral, [`certify_pclf`](@ref) decides the question exactly.

`rng` is the caller's, and `nothing` uses the default global one; pass a seeded generator when a
run has to be reproducible.

See also [`certify_pclf`](@ref).
"""
function check_pclf(pclf::PCLF, f::HybridSystems.HybridSystem; kwargs...)
    return check_pclf(pclf, UT.mode_matrices(f); kwargs...)
end

function check_pclf(
    pclf::PCLF,
    A::Vector{<:AbstractMatrix};
    nsample::Int = 20_000,
    rng = nothing,
    rtol::Float64 = 1e-6,
)
    isempty(pclf.pieces) &&
        error("The certificate carries no pieces; there is nothing to check.")
    n = size(A[1], 2)
    by_edge = Dict{Any, Float64}()

    for (u, v, label) in pclf.graph.edges
        pu, pv = _edge_pieces(pclf, u, v)
        d = piece_degree(pu)
        Am = A[Int(label)]
        worst = 0.0
        for _ in 1:nsample
            x = rng === nothing ? randn(n) : randn(rng, n)
            Vu = piece_value(pu, x)
            # A direction on which the source value vanishes carries no information about the
            # ratio, and dividing by it would manufacture one.
            Vu > 0 || continue
            worst = max(worst, piece_value(pv, Am * x) / Vu)
        end
        by_edge[(u, v, label)] = worst^(1 / d)
    end

    return _rate_report(by_edge, pclf.JSRapprox, rtol)
end

"""
    certify_pclf(pclf, f, optimizer; rtol = 1e-6)
    certify_pclf(pclf, A, optimizer; rtol = 1e-6)

Decide the contraction inequality of `pclf` exactly, for polyhedral pieces.

For a polyhedral piece the value is a gauge, so by homogeneity the rate along an edge `(s, s′, m)`
is

    sup { V_{s′}(A_m x) : V_s(x) ≤ 1 },

a maximum of linear functionals over a polytope, hence a linear program per row of `G_{s′}` and per
sign. `optimizer` is any LP-capable `JuMP` optimizer.

The return value has the same fields as [`check_pclf`](@ref), but the rates are suprema rather than
sample maxima: an empty `violated` here establishes the certificate instead of merely failing to
refute it.
"""
function certify_pclf(pclf::PCLF, f::HybridSystems.HybridSystem, optimizer; kwargs...)
    return certify_pclf(pclf, UT.mode_matrices(f), optimizer; kwargs...)
end

function certify_pclf(
    pclf::PCLF,
    A::Vector{<:AbstractMatrix},
    optimizer;
    rtol::Float64 = 1e-6,
)
    isempty(pclf.pieces) &&
        error("The certificate carries no pieces; there is nothing to certify.")
    all(p isa PolyhedralPiece for p in values(pclf.pieces)) || error(
        "certify_pclf is exact only for polyhedral pieces; use check_pclf for the others.",
    )

    n = size(A[1], 2)
    by_edge = Dict{Any, Float64}()

    for (u, v, label) in pclf.graph.edges
        pu, pv = _edge_pieces(pclf, u, v)
        M = pv.G * A[Int(label)]
        worst = 0.0
        for i in 1:size(M, 1), σ in (1.0, -1.0)
            model = JuMP.Model(optimizer)
            JuMP.set_silent(model)
            x = JuMP.@variable(model, [1:n])
            JuMP.@constraint(model, pu.G * x .<= pu.w)
            JuMP.@constraint(model, pu.G * x .>= -pu.w)
            JuMP.@objective(
                model,
                MOI.MAX_SENSE,
                σ * LinearAlgebra.dot(M[i, :], x) / pv.w[i]
            )
            JuMP.optimize!(model)
            st = JuMP.termination_status(model)
            st == MOI.OPTIMAL || error(
                "The program over edge ($u, $v, $label) terminated with status $st. An unbounded \
                 status means the sublevel set of node $u is unbounded, i.e. G is not of rank n.",
            )
            worst = max(worst, JuMP.objective_value(model))
        end
        by_edge[(u, v, label)] = worst
    end

    return _rate_report(by_edge, pclf.JSRapprox, rtol)
end

end # module
