abstract type ContinuousDomain{N, T} <: DomainType{N, T} end

"""
    ContinuousUnboundedDomain{N,T}

Struct for a basic unbounded continuous domain.
"""
struct ContinuousUnboundedDomain{N, T} <: ContinuousDomain{N, T}
    orig::SVector{N, T}
end

"""
    ContinuousBoundedDomain{N,T,B}

Struct for a basic bounded continuous domain.
"""
struct ContinuousBoundedDomain{N, T, B} <: ContinuousDomain{N, T}
    orig::SVector{N, T}
    bound::B
end

"""
    ContinuousBoundedEllipsoidDomain{N,T,S<:Grid{N,T}}

Struct for a basic bounded continuous domain formed by a finite number of ellipsoids.
"""
struct ContinuousBoundedEllipsoidDomain{N, T, B, E} <: ContinuousDomain{N, T}
    orig::SVector{N, T}
    bound::B
    ellips::Set{E}
end

function ContinuousBoundedEllipsoidDomain(orig::SVector{N, T}, bound) where {N, T}
    return ContinuousBoundedEllipsoidDomain(orig, bound, Set{UT.Ellipsoid}())
end
