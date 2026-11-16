module PathCompleteFramework

import ..Utils
const UT = Utils
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
co-complete and not complete. So the common induced on a dual graph necessarily *loses coverage*, and
the one induced on a primal graph necessarily *loses convexity*. See [`build_common_lyapunov`](@ref).
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
JSR approximation.

Dionysos does not compute certificates: it reads them as geometry. Compute one with
[PathCompleteCertificates.jl](https://github.com/dionysos-dev/PathCompleteCertificates.jl) and
convert it with `PCLF(certificate)`, which the package extension
`DionysosPathCompleteCertificatesExt` defines once both packages are loaded:

```julia
using Dionysos, PathCompleteCertificates
problem = StabilityProblem(switched_system(A))
certificate = certify(PiecewiseLinearTemplate(partitions), de_bruijn(1, 2), problem; optimizer)
pclf = Dionysos.Utils.PathCompleteFramework.PCLF(certificate)
```

A certificate can also be assembled by hand from a [`LabDigraph`](@ref) and its pieces.
"""
mutable struct PCLF
    graph::LabDigraph
    pieces::Dict{Any, AbstractPiece}
    JSRapprox::Float64
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

end # module
