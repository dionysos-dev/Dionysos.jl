module DionysosPathCompleteCertificatesExt

# Certificates from PathCompleteCertificates.jl, read as Dionysos geometry. PathCompleteCertificates
# computes the certificate; Dionysos only needs its graph, its rate, and each node function as a set
# it can cut, intersect and take pre-images of. The node functions are SetProg sets, converted here
# into the pieces `PathCompleteFramework` already knows.

import Dionysos
import PathCompleteCertificates as PCC
import LinearAlgebra

const PCF = Dionysos.Utils.PathCompleteFramework
const SetProg = PCC.SetProg

# A convex piecewise linear function, linear on each cone of a partition.
const _LinearPolytope = SetProg.Sets.Piecewise{<:Any, <:SetProg.Sets.PolarPoint}

"""
    PCLF(certificate::PathCompleteCertificates.StabilityCertificate; labels = nothing)

The certificate computed by PathCompleteCertificates.jl, as a Dionysos [`PCLF`](@ref): its graph,
its rate as `JSRapprox`, and one piece per node.

- `DualPiecewiseLinearTemplate` gives a `PolyhedralPiece(G, w)` per node, `V(x) = max_k
  |(G x)_k| / w_k`, with one row of `G` per opposite pair of the partition's vertices.
- `PiecewiseLinearTemplate` gives a `PolyhedralPiece(P, ones(m))`, `max_k |(P x)_k|`, with
  the linear pieces as rows. A Dionysos polyhedral piece is symmetric and the template's
  functions need not be, so this is `max(V(x), V(-x))`: `V` itself when it is symmetric, and
  otherwise still a certificate at the same rate, since the dynamics are linear. Its sublevel
  sets then have up to twice as many facets, and the quotient's cells follow: the dual
  template, symmetric by construction, keeps them simpler.
- `QuadraticTemplate` gives an `EllipsoidalPiece(P)`, `V(x) = x' P x`.

Nodes are numbered as in the certificate's graph, unless `labels[k]` names node `k`. The other
templates have no polyhedral or ellipsoidal geometry and are refused, as are infeasible
certificates and graphs whose edges read words: expand those with `expanded_form`.

Loading both packages also lets `MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"),
certificate)` take the certificate directly.
"""
function PCF.PCLF(certificate::PCC.StabilityCertificate; labels = nothing)
    PCC.is_feasible(certificate) || throw(
        ArgumentError("the certificate is infeasible: there is no function to convert"),
    )
    graph = PCC.get_graph(certificate)
    nodes = collect(PCC.nodes(graph))
    label(k) = labels === nothing ? k : labels[k]
    edges = [
        (label(PCC.source(e)), label(PCC.dest(e)), _mode(PCC.label(graph, e))) for
        e in PCC.edges(graph)
    ]
    U = typeof(label(first(nodes)))
    lab_digraph = PCF.LabDigraph{Int, U}(edges, Set{U}(label.(nodes)))
    pieces = Dict{Any, PCF.AbstractPiece}(
        label(k) => _piece(V) for (k, V) in zip(nodes, PCC.functions(certificate))
    )
    return PCF.PCLF(lab_digraph, pieces, Float64(certificate.rate))
end

Base.convert(::Type{PCF.PCLF}, certificate::PCC.StabilityCertificate) =
    PCF.PCLF(certificate)

_mode(mode::Integer) = Int(mode)
function _mode(word::AbstractVector)
    length(word) == 1 && return Int(only(word))
    return throw(
        ArgumentError(
            "an edge reads the word $word, and a Dionysos graph reads one mode per edge: " *
            "convert the certificate of `expanded_form(graph)` instead",
        ),
    )
end

_piece(V::SetProg.Sets.Ellipsoid) = PCF.EllipsoidalPiece(Matrix{Float64}(V.Q))

# V(x) = max_i a_i'x, the maximum of its linear pieces, and max_i |a_i'x| = max(V(x), V(-x)). The
# two are equal when V is symmetric. When it is not, the second is still a certificate at the same
# rate: the dynamics are linear, so V_d(A x) ≤ γ V_s(x) for every x gives it at -x too. One row per
# direction up to sign: a repeated row would add redundant halfspaces to every sublevel set.
function _piece(V::_LinearPolytope)
    rows = _up_to_sign([Vector{Float64}(piece.a) for piece in V.sets])
    return PCF.PolyhedralPiece(permutedims(reduce(hcat, rows)), ones(length(rows)))
end

# The polar of a piecewise linear set: V(x) = max over the partition's rays r of r'x / (a'r), with a
# the piece of a cone containing r. With r normalized, that is max_k |g_k'x| / w_k, g_k = r and
# w_k = a'r, over one ray of each opposite pair. The template imposes symmetry, so the weights of r
# and -r agree, and keeping the smaller one would give max(V(x), V(-x)) were they not to.
function _piece(V::SetProg.Sets.PolarOf{<:_LinearPolytope})
    inner = SetProg.polar(V)
    rays = Vector{Float64}[]
    weights = Float64[]
    for (piece, cone) in zip(inner.sets, inner.pieces), ray in SetProg.rays(cone)
        r = Vector{Float64}(SetProg.Polyhedra.coord(ray))
        r ./= LinearAlgebra.norm(r)
        w = LinearAlgebra.dot(piece.a, r)
        k = findfirst(q -> isapprox(q, r) || isapprox(q, -r), rays)
        if k === nothing
            push!(rays, r)
            push!(weights, w)
        else
            # A ray bounds two cones, and continuity gives both the same value up to rounding.
            # The smaller weight is the larger term of the maximum.
            weights[k] = min(weights[k], w)
        end
    end
    return PCF.PolyhedralPiece(permutedims(reduce(hcat, rays)), weights)
end

function _piece(V)
    return throw(
        ArgumentError(
            "Dionysos reads polyhedral and ellipsoidal node functions only, not $(typeof(V)): " *
            "use QuadraticTemplate, PiecewiseLinearTemplate or DualPiecewiseLinearTemplate",
        ),
    )
end

function _up_to_sign(rows)
    kept = Vector{Float64}[]
    for a in rows
        same(q) = isapprox(q, a; rtol = 1e-7) || isapprox(q, -a; rtol = 1e-7)
        any(same, kept) || push!(kept, a)
    end
    return kept
end

end # module
