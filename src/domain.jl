abstract type Domain{N,T} end

# Without S, add_set! and remove_set! where not type-stable...
struct DomainList{N,T,S<:Grid{N,T}} <: Domain{N,T}
    grid::S
    elems::Set{NTuple{N,Int}}
end

function DomainList(grid::S) where {N,S<:Grid{N}}
    return DomainList(grid, Set{NTuple{N,Int}}())
end

function add_pos!(domain::DomainList, pos)
    push!(domain.elems, pos)
end

function add_coord!(domain, x)
    add_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function add_set!(domain, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        add_pos!(domain, pos)
    end
end

function add_subset!(domain1, domain2, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain1.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain2)
        for pos in pos_iter
            if pos ∈ domain2
                add_pos!(domain1, pos)
            end
        end
    else
        for pos in enum_pos(domain2)
            if pos ∈ rectI
                add_pos!(domain1, pos)
            end
        end
    end
end

function remove_pos!(domain::DomainList, pos)
    delete!(domain.elems, pos)
end

function remove_coord!(domain, x)
    remove_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function remove_set!(domain, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain)
        for pos in pos_iter
            remove_pos!(domain, pos)
        end
    else
        for pos in enum_pos(domain)
            if pos ∈ rectI
                remove_pos!(domain, pos)
            end
        end
    end
end

function Base.union!(domain1::DomainList, domain2::DomainList)
    union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::DomainList, domain2::DomainList)
    setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::DomainList)
    empty!(domain.elems)
end

function Base.in(pos, domain::DomainList)
    return in(pos, domain.elems)
end

function Base.isempty(domain::DomainList)
    return isempty(domain.elems)
end

function Base.issubset(domain1::DomainList, domain2::DomainList)
    return issubset(domain1.elems, domain2.elems)
end

function get_ncells(domain::DomainList)
    return length(domain.elems)
end

function get_somepos(domain::DomainList)
    return first(domain.elems)
end

function enum_pos(domain::DomainList)
    return domain.elems
end
