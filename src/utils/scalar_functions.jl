abstract type ScalarFunction end

struct ZeroFunction <: ScalarFunction end

struct ConstantFunction{T} <: ScalarFunction
    value::T
end
function_value(f::ConstantFunction, x) = f.value
function Base.:+(f::ConstantFunction, g::ConstantFunction)
    return ConstantFunction(f.value + g.value)
end

struct AffineFunction{T}
    a::Vector{T}
    β::T
end

function function_value(f::AffineFunction, x)
    return sum(f.a .* x) + f.β
end
function Base.isapprox(f::AffineFunction, g::AffineFunction; kws...)
    return isapprox(f.a, g.a; kws...) && isapprox(f.β, g.β; kws...)
end

struct QuadraticControlFunction{T, MT <: AbstractMatrix{T}} <: ScalarFunction
    Q::MT
end
function function_value(f::QuadraticControlFunction, x)
    return x'f.Q * x
end

struct PolyhedralFunction{T} <: ScalarFunction
    lower_bound::T
    pieces::Vector{AffineFunction{T}}
    domain::Polyhedra.Intersection{T, Vector{T}, Int}
end
_inf(T::Type{<:AbstractFloat}) = typemax(T)
_inf(T::Type) = error("No infinite value for type $T")
function function_value(f::PolyhedralFunction{T}, x) where {T}
    if !(x in f.domain)
        return _inf(T)
    end
    return mapreduce(piece -> function_value(piece, x), max, f.pieces; init = f.lower_bound)
end

function Base.:+(c::ConstantFunction, p::PolyhedralFunction)
    return PolyhedralFunction(c.value + p.lower_bound, p.pieces, p.domain)
end

Base.:+(::ZeroFunction, f::Union{ConstantFunction, PolyhedralFunction}) = f

abstract type ScalarControlFunction end
struct ConstantControlFunction{T} <: ScalarControlFunction
    value::T
end
function_value(f::ConstantControlFunction, x, u) = f.value

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
function function_value(f::QuadraticStateControlFunction, x, u)
    return x'f.Q * x + u'f.R * u + 2 * (x'f.N * u + x'f.q + u'f.r) + f.v
end
function get_full_psd_matrix(f::QuadraticStateControlFunction)
    return [
        f.Q f.N f.q
        f.N' f.R f.r
        f.q' f.r' f.v
    ]
end
