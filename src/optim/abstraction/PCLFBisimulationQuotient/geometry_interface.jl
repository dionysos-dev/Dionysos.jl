# ============================================================
# Utilities
# ============================================================

const Poly = HPolytope

# Try to coerce a LazySets polyhedral object to an `HPolytope`.
function _as_hpolytope(P)
    P isa HPolytope && return P
    try
        return HPolytope(constraints_list(P))
    catch
        error("Cannot convert object of type $(typeof(P)) to HPolytope.")
    end
end

function set_intersection(P::Poly, Q::Poly)
    return HPolytope(vcat(constraints_list(P), constraints_list(Q)))
end

function is_nonempty_set(P::Poly)
    return !isempty(P)
end

"""
    preimage_linear(P, A)

Return the linear preimage `{x : A*x ∈ P}`.
If `P = {y : H y ≤ h}`, then `Pre(P,A) = {x : H A x ≤ h}`.
"""
function preimage_linear(P::Poly, A::AbstractMatrix)
    new_cons = HalfSpace[]
    for c in constraints_list(P)
        push!(new_cons, HalfSpace(A' * c.a, c.b))
    end
    return HPolytope(new_cons)
end

# ============================================================
# Polyhedral difference decomposition
# ============================================================

"""
    set_difference_decompose(P, Q; atol=0.0)

Decompose `P \\ Q` into a finite union of H-polytopes using a sequential
halfspace decomposition.

If `Q = ⋂_{j=1}^m H_j`, then
`P \\ Q = ⋃_{j=1}^m ( P ∩ H_1 ∩ ... ∩ H_{j-1} ∩ H_j^c )`.

Here `H_j^c` is represented by the closed opposite halfspace
`-a'x <= -b - atol` for `H_j = {x : a'x <= b}`.

Notes:
- `atol = 0.0` keeps the boundary and may produce overlaps on boundaries.
- `atol > 0` separates boundaries slightly but changes the set by a small
  tolerance.
- This is the main place you may later replace by an exact Baotić-style
  decomposition.
"""
function set_difference_decompose(P::Poly, Q::Poly; atol::Float64 = 0.0)
    qcons = constraints_list(Q)
    pcons = constraints_list(P)

    pieces = Poly[]
    prefix = HalfSpace[]

    for c in qcons
        # Complement of a'x <= b is approximated by a'x >= b + atol,
        # i.e. -a'x <= -(b + atol)
        comp = HalfSpace(-c.a, -(c.b + atol))
        piece = HPolytope(vcat(pcons, prefix, [comp]))
        if is_nonempty_set(piece)
            push!(pieces, piece)
        end
        push!(prefix, c)
    end

    return pieces
end

"""
    set_difference_decompose(P, Qs::Vector{Poly}; atol=0.0)

Iteratively decompose `P \\ (⋃ Qs)` as a union of polytopes.
"""
function set_difference_decompose(P::Poly, Qs::AbstractVector{<:Poly}; atol::Float64 = 0.0)
    parts = Poly[P]
    for Q in Qs
        new_parts = Poly[]
        for R in parts
            append!(new_parts, set_difference_decompose(R, Q; atol = atol))
        end
        parts = new_parts
    end
    return parts
end
