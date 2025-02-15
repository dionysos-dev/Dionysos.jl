"""
    CustomList{N, T} <: DomainType{N, T}

A flexible, generic domain representation that stores elements as a list of `SVector{N, T}`.  
Useful for managing discrete sets of points in an `N`-dimensional space.
"""
struct CustomList{N, T} <: DomainType{N, T}
    elems::Vector{SVector{N, T}}
end

enum_elems(domain::CustomList) = domain.elems
get_ncells(domain::CustomList) = length(domain.elems)
get_elem_by_index(domain::CustomList, index) = domain.elems[index]
Base.isempty(domain::CustomList) = isempty(domain.elems)
Base.union!(domain1::CustomList, domain2::CustomList) = union!(domain1.elems, domain2.elems)
Base.setdiff!(domain1::CustomList, domain2::CustomList) =
    setdiff!(domain1.elems, domain2.elems)
Base.empty!(domain::CustomList) = empty!(domain.elems)

convert_to_custom_domain(domain::CustomList) = domain
