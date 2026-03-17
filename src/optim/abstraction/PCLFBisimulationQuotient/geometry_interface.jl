# ============================================================
# Polytopes
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

function clean_poly(P::Poly)
    return remove_redundant_constraints(P)
end

function set_intersection(P::Poly, Q::Poly)
    return clean_poly(HPolytope(vcat(constraints_list(P), constraints_list(Q))))
    # intersection = HPolytope(vcat(constraints_list(P), constraints_list(Q)))
    # if !isempty(intersection)
    #     return clean_poly(intersection)
    # else
    #     return EmptySet{Float64}(dim(P))
    # end
end

is_nonempty_set(::LazySets.EmptySet) = false
is_nonempty_set(P) = !isempty(P)

# ============================================================
# Semi-linear sets = finite unions of polytopes
# ============================================================

mutable struct SemiLinearSet
    parts::Vector{Poly}

    function SemiLinearSet(parts::Vector{Poly})
        return new(parts)
    end
end

SemiLinearSet() = SemiLinearSet(Poly[])

function SemiLinearSet(p::Poly)
    return SemiLinearSet(Poly[_as_hpolytope(p)])
end

function SemiLinearSet(parts::AbstractVector)
    out = Poly[]
    for P in parts
        push!(out, _as_hpolytope(P))
    end
    return SemiLinearSet(out)
end

Base.length(S::SemiLinearSet) = length(S.parts)
Base.isempty(S::SemiLinearSet) = isempty(S.parts)
Base.iterate(S::SemiLinearSet, st...) = iterate(S.parts, st...)

function Base.show(io::IO, S::SemiLinearSet)
    return print(io, "SemiLinearSet($(length(S.parts)) parts)")
end

function normalize_semilinear(S::SemiLinearSet)
    keep = Poly[]
    for Pk in S.parts
        if is_nonempty_set(Pk)
            push!(keep, _as_hpolytope(Pk))
        end
    end
    return SemiLinearSet(keep)
end

# ============================================================
# Semi-linear sets : Nonemptiness
# ============================================================

function is_nonempty_set(S::SemiLinearSet)
    return any(is_nonempty_set(P) for P in S.parts)
end

# ============================================================
# Semi-linear sets : Intersection
# ============================================================

function set_intersection(S1::SemiLinearSet, S2::SemiLinearSet)
    out = Poly[]
    for P1 in S1.parts, P2 in S2.parts
        I = set_intersection(P1, P2)
        if is_nonempty_set(I)
            push!(out, _as_hpolytope(I))
        end
    end
    return SemiLinearSet(out)
end

function set_intersection(S::SemiLinearSet, P0::Poly)
    out = Poly[]
    for P1 in S.parts
        I = set_intersection(P1, P0)
        if is_nonempty_set(I)
            push!(out, _as_hpolytope(I))
        end
    end
    return SemiLinearSet(out)
end

function set_intersection(P0::Poly, S::SemiLinearSet)
    return set_intersection(S, P0)
end

# ============================================================
# Semi-linear sets : Difference
# ============================================================

function set_difference_decompose(S1::SemiLinearSet, S2::SemiLinearSet; atol::Float64 = 0.0)
    current = copy(S1.parts)
    for Q in S2.parts
        new_current = Poly[]
        for P1 in current
            append!(new_current, set_difference_decompose(P1, Q; atol = atol))
        end
        current = new_current
    end
    return SemiLinearSet(current)
end

function set_difference_decompose(S::SemiLinearSet, P0::Poly; atol::Float64 = 0.0)
    out = Poly[]
    for P1 in S.parts
        append!(out, set_difference_decompose(P1, P0; atol = atol))
    end
    return SemiLinearSet(out)
end

function set_difference_decompose(P::Poly, Q::Poly; atol::Float64 = 0.0)
    qcons = constraints_list(Q)
    pcons = constraints_list(P)

    pieces = Poly[]
    prefix = HalfSpace[]

    for c in qcons
        comp = HalfSpace(-c.a, -(c.b + atol))
        piece = clean_poly(HPolytope(vcat(pcons, prefix, [comp])))
        if is_nonempty_set(piece)
            push!(pieces, _as_hpolytope(piece))
        end
        push!(prefix, c)
    end
    return pieces
end

# ============================================================
# Semi-linear sets : Preimage under linear map
# ============================================================

function preimage_linear(S::SemiLinearSet, A::AbstractMatrix)
    out = Poly[]
    for P1 in S.parts
        Ppre = preimage_linear(P1, A)
        if is_nonempty_set(Ppre)
            push!(out, _as_hpolytope(Ppre))
        end
    end
    return SemiLinearSet(out)
end

function preimage_linear(P::Poly, A::AbstractMatrix)
    new_cons = HalfSpace[]
    for c in constraints_list(P)
        push!(new_cons, HalfSpace(A' * c.a, c.b))
    end
    return clean_poly(HPolytope(new_cons))
end

function preimage_linear_parts(S::SemiLinearSet, A::AbstractMatrix)
    parts = Poly[]
    for P in S
        preP = preimage_linear(P, A)
        if is_nonempty_set(preP)
            push!(parts, preP)
        end
    end
    return parts
end

@recipe function f(
    S::SemiLinearSet;
    fillcolor = :blue,
    linecolor = :blue,
    fillalpha = 0.25,
    linealpha = 1.0,
    linewidth = 1.0,
    show_label = false,
)
    for (k, P) in enumerate(S.parts)
        @series begin
            fillcolor := fillcolor
            linecolor := linecolor
            fillalpha := fillalpha
            linealpha := linealpha
            linewidth := linewidth
            label := (show_label && k == 1) ? "SemiLinearSet" : ""
            P
        end
    end
end
