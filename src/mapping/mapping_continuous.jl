"""
    MappingContinuousGrid(C::DO.ContinuousDomain, D::DomainList{N,T,S<:DO.Grid{N,T}}

Mapping containing methods to convert elements in a `ContinuousDomain` to a grid-based `DomainList`
and vice-versa.
"""
struct MappingContinuousGrid{C<:DO.ContinuousDomain{N,T}, D<:DO.DomainList{N,T,S<:DO.Grid{N,T}}}
    continuousdomain::C
    griddomain::D
    cont2grid::Function
    grid2cont::Function
end

function MappingContinuousGrid(c<:DO.ContinuousDomain{N,T}, g<:DO.DomainList{N,T,S<:DO.Grid{N,T}}) where {N,S<:Grid{N}}
    if (g.grid <: DO.GridEllipsoidalRectangular)
        cont2grid(x) = DO.get_ellip(g.grid, x)
        grid2cont(xi) = xi.c
    else
        cont2grid(x) = DO.get_rec(g.grid, x)
        grid2cont(xi) = UT.get_center(xi)
    end
    return MappingContinuousGrid(c, g, cont2grid, grid2cont)
end



"""
    MappingContinuousEllipsoid{C<:DO.ContinuousDomain{N,T}, CE<:DO.ContinuousBoundedEllipsoidDomain{N,T,B,E}}

Mapping containing methods to convert elements in a `ContinuousDomain` to a `ContinuousBoundedEllipsoidDomain`
and vice-versa.
"""
struct MappingContinuousEllipsoid{C<:DO.ContinuousDomain{N,T}, CE<:DO.ContinuousBoundedEllipsoidDomain{N,T,B,E}}
    continuousdomain::C
    griddomain::D
    cont2elli::Function
    elli2cont::Function
end


function MappingContinuousEllipsoid(c<:DO.ContinuousDomain{N,T}, e<:DO.ContinuousBoundedEllipsoidDomain{N,T,B,E}) 
    cont2elli(x) = e.ellips[findfirst(x .∈ e.ellips)] 
    elli2cont(xi) = xi.c
    return MappingContinuousGrid(c, g, cont2grid, grid2cont)
end
