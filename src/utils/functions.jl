# Structured scalar functions: plain-data (hence serializable, e.g. to JLD2)
# representations of costs and value functions, with algebra (`+`) and
# solver-side dispatch on the concrete type.
#
# The evaluation contract is *callable*: a state cost is any callable `f(x)`,
# a transition cost any callable `f(x, u)` — so a plain `Function` or closure
# is always a valid cost; the structs below are the representations that also
# survive serialization and structural introspection.

"Root of state functions `f(x)`."
abstract type ScalarFunction end

"Root of state-input functions `f(x, u)`."
abstract type ScalarControlFunction end

struct ZeroFunction <: ScalarFunction end
(f::ZeroFunction)(x) = 0.0
(f::ZeroFunction)(x, u) = 0.0

struct ConstantFunction{T} <: ScalarFunction
    value::T
end
(f::ConstantFunction)(x) = f.value
function Base.:+(f::ConstantFunction, g::ConstantFunction)
    return ConstantFunction(f.value + g.value)
end

struct AffineFunction{T} <: ScalarFunction
    a::Vector{T}
    β::T
end
(f::AffineFunction)(x) = sum(f.a .* x) + f.β
function Base.isapprox(f::AffineFunction, g::AffineFunction; kws...)
    return isapprox(f.a, g.a; kws...) && isapprox(f.β, g.β; kws...)
end

"Quadratic form `x ↦ x'Qx`."
struct QuadraticFunction{T, MT <: AbstractMatrix{T}} <: ScalarFunction
    Q::MT
end
(f::QuadraticFunction)(x) = x' * f.Q * x

"Piecewise-affine convex function `x ↦ max(lower_bound, maxᵢ pieceᵢ(x))` on `domain` (`∞` outside)."
struct PolyhedralFunction{T} <: ScalarFunction
    lower_bound::T
    pieces::Vector{AffineFunction{T}}
    domain::Polyhedra.Intersection{T, Vector{T}, Int}
end
_inf(T::Type{<:AbstractFloat}) = typemax(T)
_inf(T::Type) = error("No infinite value for type $T")
function (f::PolyhedralFunction{T})(x) where {T}
    if !(x in f.domain)
        return _inf(T)
    end
    return mapreduce(piece -> piece(x), max, f.pieces; init = f.lower_bound)
end

function Base.:+(c::ConstantFunction, p::PolyhedralFunction)
    return PolyhedralFunction(c.value + p.lower_bound, p.pieces, p.domain)
end

Base.:+(::ZeroFunction, f::Union{ConstantFunction, PolyhedralFunction}) = f

"""
    BlackBoxFunction{F}

Arbitrary user-supplied state function `f(x)`, in the spirit of
MathematicalSystems' `BlackBox` maps: full modeling freedom, but unlike the
structured functions it does not serialize reliably (closures, e.g. to JLD2).
"""
struct BlackBoxFunction{F} <: ScalarFunction
    f::F
end
(g::BlackBoxFunction)(x) = g.f(x)

struct ConstantControlFunction{T} <: ScalarControlFunction
    value::T
end
(f::ConstantControlFunction)(x, u) = f.value

"""
    BlackBoxControlFunction{F}

Arbitrary user-supplied state-input function `f(x, u)`; same freedom and
serialization caveat as [`BlackBoxFunction`](@ref).
"""
struct BlackBoxControlFunction{F} <: ScalarControlFunction
    f::F
end
(g::BlackBoxControlFunction)(x, u) = g.f(x, u)

"""
    QuadraticStateControlFunction{T, MT<:AbstractMatrix{T}}

Quadratic function on state and input defined as
`x'Qx + u'Ru + 2x'Nu + 2x'q + 2u'r + v`
"""
struct QuadraticStateControlFunction{T, MT <: AbstractMatrix{T}, AT <: AbstractArray{T}} <:
       ScalarControlFunction
    Q::MT
    R::MT
    N::MT
    q::AT
    r::AT
    v::T
end
function (f::QuadraticStateControlFunction)(x, u)
    return x' * f.Q * x + u' * f.R * u + 2 * (x' * f.N * u + x' * f.q + u' * f.r) + f.v
end
function get_full_psd_matrix(f::QuadraticStateControlFunction)
    return [
        f.Q f.N f.q
        f.N' f.R f.r
        f.q' f.r' f.v
    ]
end

# Canonical `[x; u; 1]` PSD cost matrix, so every SDP/MIQP consumer converts a
# transition cost the same way instead of re-deriving the dispatch. A raw matrix
# is assumed to already be in this form.
_cost_matrix(M::AbstractMatrix) = M
_cost_matrix(f::QuadraticStateControlFunction) = get_full_psd_matrix(f)
