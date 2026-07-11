# ============================================================
# Semi-linear sets = finite unions of polytopes, represented natively as a
# `LazySets.UnionSetArray` of `HPolytope`s. LazySets provides membership,
# emptiness, dimension, concrete `intersection`/`isdisjoint` and plotting.
# Dionysos keeps only what LazySets lacks:
#   - `set_difference_decompose`, the tolerance-based difference kernel: cut
#     halfspaces are shifted by `atol`, so pieces stay disjoint and degenerate
#     slivers are dropped — piece counts (and hence the PCLF-bisimulation
#     runtime) depend on this, which is why we do not use the exact
#     `LazySets.difference`;
#   - `disjointify` and the union volume built on it (LazySets' union volume
#     uses inclusion–exclusion, a different algorithm);
#   - `preimage_linear`, the halfspace preimage under a linear map.
# ============================================================

function _as_hpolytope(P)
    P isa LazySets.HPolytope && return P
    try
        return LazySets.HPolytope(LazySets.constraints_list(P))
    catch
        error("Cannot convert object of type $(typeof(P)) to HPolytope.")
    end
end

function clean_poly(P::LazySets.HPolytope)
    try
        return LazySets.remove_redundant_constraints(P)
    catch
        @warn "clean_poly: remove_redundant_constraints failed; returning input polytope" maxlog =
            1
        return P
    end
end

"""
    SemiLinearSet{T, PT}

A finite union of polytopes: alias of `LazySets.UnionSetArray` restricted to
`HPolytope` parts, so the whole LazySets API (`∈`, `isempty`, `dim`,
`intersection`, plotting) applies. Build with [`semilinear_set`](@ref); the
parts live in `S.array`.
"""
const SemiLinearSet{T, PT <: LazySets.HPolytope{T}} = LazySets.UnionSetArray{T, PT}

"""
    semilinear_set(parts) -> SemiLinearSet

The union of the given polyhedral sets, each converted to `HPolytope`.
"""
function semilinear_set(parts::AbstractVector)
    out = [_as_hpolytope(P) for P in parts]
    isempty(out) && return semilinear_set()
    ET = eltype(out)
    isconcretetype(ET) && return LazySets.UnionSetArray(out)
    T = _num_type(first(out))
    return LazySets.UnionSetArray(convert(Vector{LazySets.HPolytope{T}}, out))
end
semilinear_set() = LazySets.UnionSetArray(LazySets.HPolytope{Float64, Vector{Float64}}[])
semilinear_set(P::LazySets.LazySet) = semilinear_set([P])

# Drop empty parts (each check is one LP feasibility problem).
function normalize_semilinear(S::SemiLinearSet)
    return semilinear_set([P for P in S.array if !isempty(P)])
end

# Concrete polytope intersection, always an `HPolytope` or `EmptySet`. The
# fast path is `LazySets.intersection` (in 2-D an LP-free polygon algorithm —
# refinement loops depend on that speed), but on degenerate lower-dimensional
# intersections it can error internally (e.g. `element(::Line2D)`); the
# fallback then stacks constraints and prunes with one LP per constraint.
function poly_intersection(P::LazySets.HPolytope, Q::LazySets.HPolytope)
    I = try
        LazySets.intersection(P, Q)
    catch
        nothing
    end
    if I isa LazySets.HPolytope || I isa LazySets.EmptySet
        return I
    end
    raw =
        LazySets.HPolytope(vcat(LazySets.constraints_list(P), LazySets.constraints_list(Q)))
    isempty(raw) && return LazySets.EmptySet{_num_type(P)}(LazySets.dim(P))
    return clean_poly(raw)
end

# Nonempty parts of `S ∩ P0`, as a vector of `HPolytope`s (`poly_intersection`
# already prunes empty results, so no extra feasibility LP is spent here).
function poly_intersection_parts(S::SemiLinearSet, P0::LazySets.HPolytope)
    out = LazySets.HPolytope[]
    for P1 in S.array
        I = poly_intersection(P1, P0)
        I isa LazySets.EmptySet || push!(out, I)
    end
    return out
end

# ============================================================
# Volume (via a Polyhedra backend; overlapping parts are disjointified first)
# ============================================================

function disjointify(S::SemiLinearSet; atol::Float64 = 0.0)
    pieces = LazySets.HPolytope[]

    for P in S.array
        current = [_as_hpolytope(P)]

        for Q in pieces
            new_current = LazySets.HPolytope[]
            for C in current
                append!(new_current, set_difference_decompose(C, Q; atol = atol))
            end
            current = new_current
            isempty(current) && break
        end

        append!(pieces, current)
    end

    return semilinear_set(pieces)
end

function get_volume(
    S::SemiLinearSet;
    backend = nothing,
    assume_disjoint::Bool = false,
    atol::Float64 = 0.0,
)
    Sd = assume_disjoint ? S : disjointify(S; atol = atol)
    return sum(LazySets.volume(P; backend = backend) for P in Sd.array)
end

# ============================================================
# Difference — tolerance-based carve-out decomposition (see header)
# ============================================================

function set_difference_decompose(
    S1::SemiLinearSet,
    S2::SemiLinearSet;
    atol::Float64 = 1e-6,
)
    current = copy(S1.array)
    for Q in S2.array
        new_current = LazySets.HPolytope[]
        for P1 in current
            append!(new_current, set_difference_decompose(P1, Q; atol = atol))
        end
        current = new_current
    end
    return current
end

function set_difference_decompose(
    S::SemiLinearSet,
    P0::LazySets.HPolytope;
    atol::Float64 = 1e-6,
)
    out = LazySets.HPolytope[]
    for P1 in S.array
        append!(out, set_difference_decompose(P1, P0; atol = atol))
    end
    return out
end

function set_difference_decompose(
    P::LazySets.HPolytope,
    S::SemiLinearSet;
    atol::Float64 = 1e-6,
)
    current = [P]

    for Q in S.array
        new_current = LazySets.HPolytope[]
        for P1 in current
            append!(new_current, set_difference_decompose(P1, Q; atol = atol))
        end
        current = new_current
    end

    return current
end

function set_difference_decompose(
    P::LazySets.HPolytope,
    Q::LazySets.HPolytope;
    atol::Float64 = 1e-6,
)
    qcons = LazySets.constraints_list(Q)
    pcons = LazySets.constraints_list(P)

    pieces = LazySets.HPolytope[]
    prefix = LazySets.HalfSpace[]

    for c in qcons
        comp = LazySets.HalfSpace(-c.a, -(c.b + atol))
        piece = clean_poly(LazySets.HPolytope(vcat(pcons, prefix, [comp])))
        if !isempty(piece)
            push!(pieces, piece)
        end
        push!(prefix, c)
    end
    return pieces
end

# ============================================================
# Preimage under a linear map: {x : Ax ∈ P} (no LazySets equivalent)
# ============================================================

function preimage_linear(S::SemiLinearSet, A::AbstractMatrix)
    out = LazySets.HPolytope[]
    for P1 in S.array
        Ppre = preimage_linear(P1, A)
        if !isempty(Ppre)
            push!(out, Ppre)
        end
    end
    return semilinear_set(out)
end

function preimage_linear(P::LazySets.HPolytope, A::AbstractMatrix)
    new_cons = LazySets.HalfSpace[]
    for c in LazySets.constraints_list(P)
        push!(new_cons, LazySets.HalfSpace(A' * c.a, c.b))
    end
    return clean_poly(LazySets.HPolytope(new_cons))
end

function preimage_linear_parts(S::SemiLinearSet, A::AbstractMatrix)
    parts = LazySets.HPolytope[]
    for P in S.array
        preP = preimage_linear(P, A)
        if !isempty(preP)
            push!(parts, preP)
        end
    end
    return parts
end
