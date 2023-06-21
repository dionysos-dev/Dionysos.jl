"""
    MappingContinuousGrid(C::DO.ContinuousDomain, D::DomainList{N,T,S<:DO.Grid{N,T}}

Mapping containing methods to convert elements in a `ContinuousDomain` to a grid-based `DomainList`
and vice-versa.
"""
struct MappingContinuousGrid{
    N,
    T,
    F1 <: Function,
    F2 <: Function,
    S <: DO.Grid{N, T},
    C <: DO.ContinuousDomain{N, T},
    D <: DO.DomainList{N, T, S},
}
    continuousdomain::C
    griddomain::D
    cont2grid::F1
    grid2cont::F2
end

function MappingContinuousGrid(
    c::C,
    g::DO.DomainList{N, T, S},
) where {N, T, S <: DO.Grid{N, T}, C <: DO.ContinuousDomain{N, T}}
    is_ellipsoid = g.grid isa DO.GridEllipsoidalRectangular

    cont2grid(x) = is_ellipsoid ? DO.get_elem_by_pos(g.grid, x) : DO.get_rec(g.grid, x)
    grid2cont(xi) = is_ellipsoid ? xi.c : UT.get_center(xi)

    return MappingContinuousGrid(c, g, cont2grid, grid2cont)
end

"""
    MappingContinuousEllipsoid{C<:DO.ContinuousDomain{N,T}, CE<:DO.ContinuousBoundedEllipsoidDomain{N,T,B,E}}

Mapping containing methods to convert elements in a `ContinuousDomain` to a `ContinuousBoundedEllipsoidDomain`
and vice-versa.
"""
struct MappingContinuousEllipsoid{
    N,
    T,
    B,
    E,
    F1 <: Function,
    F2 <: Function,
    C <: DO.ContinuousDomain{N, T},
    G <: DO.ContinuousBoundedEllipsoidDomain{N, T, B, E},
}
    continuousdomain::C
    ellidomain::G
    cont2elli::F1
    elli2cont::F2
end

function MappingContinuousEllipsoid(
    c::C,
    e::DO.ContinuousBoundedEllipsoidDomain{N, T, B, E},
) where {N, T, B, E, C <: DO.ContinuousDomain{N, T}}
    function cont2elli(x)
        idx = findfirst([x ∈ ell for ell in e.ellips])
        return idx === nothing ? idx : e.ellips[idx]
    end
    elli2cont(xi) = xi.c
    return MappingContinuousEllipsoid(c, e, cont2elli, elli2cont)
end
