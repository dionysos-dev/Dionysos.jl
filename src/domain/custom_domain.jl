
struct CustomList{N,T} <: DomainType{N,T}
    elems::Vector{SVector{N,T}}
end

# function CustomList(elems::Vector{SVector{N,T}}) where {N,T}
#     return CustomList(elems)
# end

function enum_pos(domain::CustomList)
    return domain.elems
end

function Base.isempty(domain::CustomList)
    return isempty(domain.elems)
end

function Base.union!(domain1::CustomList, domain2::CustomList)
    union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::CustomList, domain2::CustomList)
    setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::CustomList)
    empty!(domain.elems)
end

function get_ncells(domain::CustomList)
    return length(domain.elems)
end

# function get_coord(domain::CustomList, rec)
#     return U.center(rec)
# end
