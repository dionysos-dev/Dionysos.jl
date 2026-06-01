# ============================================================
# Polytopes
# ============================================================

const Poly = LazySets.HPolytope

function _as_hpolytope(P)
    P isa Poly && return P
    try
        return Poly(LazySets.constraints_list(P))
    catch
        error("Cannot convert object of type $(typeof(P)) to HPolytope.")
    end
end

function clean_poly(P::Poly)
    try
        return LazySets.remove_redundant_constraints(P)
    catch
        println("Clean failed")
        return P
    end
end

function set_intersection(P::Poly, Q::Poly)
    return clean_poly(
        Poly(vcat(LazySets.constraints_list(P), LazySets.constraints_list(Q))),
    )
end

is_nonempty_set(::LazySets.EmptySet) = false
is_nonempty_set(P) = !isempty(P)

function get_volume(P::Poly; backend = nothing)
    backend === nothing && error(
        "No polyhedral backend provided. Set `polyhedra_backend`, e.g. CDDLib.Library().",
    )
    ph = LazySets.polyhedron(P; backend = backend)
    return Polyhedra.volume(ph)
end

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
function dim(P::Poly)
    cons = LazySets.constraints_list(P)
    isempty(cons) && error("Cannot infer dimension from empty constraint list.")
    return length(first(cons).a)
end
function dim(S::SemiLinearSet)
    isempty(S) && error("Cannot infer dimension of an empty SemiLinearSet.")
    return dim(first(S.parts))
end

function Base.in(x::AbstractVector, S::SemiLinearSet)
    return any(Base.in(x, P) for P in S.parts)
end

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
# Semi-linear sets : Volume
# ============================================================

function disjointify(S::SemiLinearSet; atol::Float64 = 0.0)
    pieces = Poly[]

    for P in S.parts
        current = Poly[_as_hpolytope(P)]

        for Q in pieces
            new_current = Poly[]
            for C in current
                append!(new_current, set_difference_decompose(C, Q; atol = atol))
            end
            current = new_current
            isempty(current) && break
        end

        append!(pieces, current)
    end

    return SemiLinearSet(pieces)
end

function get_volume(
    S::SemiLinearSet;
    backend = nothing,
    assume_disjoint::Bool = false,
    atol::Float64 = 0.0,
)
    Sd = assume_disjoint ? S : disjointify(S; atol = atol)
    return sum(get_volume(P; backend = backend) for P in Sd.parts)
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

function hyperrectangle_to_hpolytope(X::LazySets.Hyperrectangle)
    c = LazySets.center(X)
    r = LazySets.radius_hyperrectangle(X)

    n = length(c)
    cons = LazySets.HalfSpace[]

    for i in 1:n
        ei = zeros(n)
        ei[i] = 1.0

        # x_i ≤ c_i + r_i
        push!(cons, LazySets.HalfSpace(ei, c[i] + r[i]))

        # -x_i ≤ -(c_i - r_i)
        push!(cons, LazySets.HalfSpace(-ei, -(c[i] - r[i])))
    end

    return LazySets.HPolytope(cons)
end

function set_intersection(P::Poly, X::LazySets.Hyperrectangle)
    return set_intersection(P, hyperrectangle_to_hpolytope(X))
end

# ============================================================
# Semi-linear sets : Difference
# ============================================================

function set_difference_decompose(
    S1::SemiLinearSet,
    S2::SemiLinearSet;
    atol::Float64 = 1e-6,
)
    current = copy(S1.parts)
    for Q in S2.parts
        new_current = Poly[]
        for P1 in current
            append!(new_current, set_difference_decompose(P1, Q; atol = atol))
        end
        current = new_current
    end
    return current
end

function set_difference_decompose(S::SemiLinearSet, P0::Poly; atol::Float64 = 1e-6)
    out = Poly[]
    for P1 in S.parts
        append!(out, set_difference_decompose(P1, P0; atol = atol))
    end
    return out
end

function set_difference_decompose(P::Poly, S::SemiLinearSet; atol::Float64 = 1e-6)
    current = [P]

    for Q in S.parts
        new_current = Poly[]
        for P1 in current
            append!(new_current, set_difference_decompose(P1, Q; atol = atol))
        end
        current = new_current
    end

    return current
end

function set_difference_decompose(P::Poly, Q::Poly; atol::Float64 = 1e-6)
    qcons = LazySets.constraints_list(Q)
    pcons = LazySets.constraints_list(P)

    pieces = Poly[]
    prefix = LazySets.HalfSpace[]

    for c in qcons
        comp = LazySets.HalfSpace(-c.a, -(c.b + atol))
        piece = clean_poly(Poly(vcat(pcons, prefix, [comp])))
        if is_nonempty_set(piece)
            push!(pieces, piece)
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
    new_cons = LazySets.HalfSpace[]
    for c in LazySets.constraints_list(P)
        push!(new_cons, LazySets.HalfSpace(A' * c.a, c.b))
    end
    return clean_poly(Poly(new_cons))
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
        try
            P = LazySets.remove_redundant_constraints(P)

            @series begin
                fillcolor := fillcolor
                linecolor := linecolor
                fillalpha := fillalpha
                linealpha := linealpha
                linewidth := linewidth
                label := (show_label && k == 1) ? "SemiLinearSet" : ""
                #               P
            end

        catch err
            @warn "Skipping problematic polytope during plotting" exception=(
                err,
                catch_backtrace(),
            )
            continue
        end
    end
end
