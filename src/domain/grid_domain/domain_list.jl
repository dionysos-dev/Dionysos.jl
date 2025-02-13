"""
    DomainList{N,T,S<:Grid{N,T}}

Struct for a basic domain based on a `Grid`.
"""
struct DomainList{N, T, S <: Grid{N, T}} <: GridDomainType{N, T}
    grid::S
    elems::Set{NTuple{N, Int}}
end

function DomainList(grid::S) where {N, S <: Grid{N}}
    return DomainList(grid, Set{NTuple{N, Int}}())
end

get_grid(domain::DomainList) = domain.grid
enum_pos(domain::DomainList) = domain.elems
get_ncells(domain::DomainList) = length(domain.elems)
get_somepos(domain::DomainList) = first(domain.elems)
Base.isempty(domain::DomainList) = isempty(domain.elems)
Base.in(pos, domain::DomainList) = in(pos, domain.elems)
Base.issubset(domain1::DomainList, domain2::DomainList) =
    issubset(domain1.elems, domain2.elems)
add_pos!(domain::DomainList, pos) = push!(domain.elems, pos)
Base.union!(domain1::DomainList, domain2::DomainList) = union!(domain1.elems, domain2.elems)
Base.setdiff!(domain1::DomainList, domain2::DomainList) =
    setdiff!(domain1.elems, domain2.elems)
Base.empty!(domain::DomainList) = empty!(domain.elems)
remove_pos!(domain::DomainList, pos) = delete!(domain.elems, pos)
