abstract type DomainType{N, T} end

function _invInclMode(incl_mode::INCL_MODE)
    return incl_mode == OUTER ? INNER : OUTER
end

# Without S, add_set! and remove_set! where not type-stable...

"""
    DomainList{N,T,S<:Grid{N,T}}

Struct for a basic domain based on a `Grid`
"""
struct DomainList{N, T, S <: Grid{N, T}} <: DomainType{N, T}
    grid::S
    elems::Set{NTuple{N, Int}}
end

"""
    DomainList(grid::S) where {N,S<:Grid{N}}

Return a new DomainList
"""
function DomainList(grid::S) where {N, S <: Grid{N}}
    return DomainList(grid, Set{NTuple{N, Int}}())
end

function get_grid(domain::DomainList)
    return domain.grid
end

function add_pos!(domain::DomainList, pos)
    return push!(domain.elems, pos)
end

function add_coord!(domain, x)
    return add_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function add_set!(domain, setMinus::UT.LazySetMinus, incl_mode::INCL_MODE)
    add_set!(domain, setMinus.A, incl_mode)
    return remove_set!(domain, setMinus.A ∩ setMinus.B, _invInclMode(incl_mode))
end

function add_set!(domain, unionSetArray::UT.LazyUnionSetArray, incl_mode::INCL_MODE)
    for set in unionSetArray.sets
        add_set!(domain, set, incl_mode)
    end
end

function add_set!(domain, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        add_pos!(domain, pos)
    end
end

function get_subset_pos(domain::DomainList, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    posL = []
    for pos in pos_iter
        if pos ∈ domain
            push!(posL, pos)
        end
    end
    return posL
end

function add_subset!(domain1, domain2, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
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
    return delete!(domain.elems, pos)
end

function remove_coord!(domain, x)
    return remove_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function remove_set!(domain, unionSetArray::UT.LazyUnionSetArray, incl_mode::INCL_MODE)
    for set in unionSetArray.sets
        remove_set!(domain, set, incl_mode)
    end
end

function remove_set!(domain, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
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
    return union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::DomainList, domain2::DomainList)
    return setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::DomainList)
    return empty!(domain.elems)
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

function crop_to_domain(domain::DomainList, list)
    return list ∩ enum_pos(domain)
end

function get_coord(domain::DomainType, pos)
    return get_coord_by_pos(domain.grid, pos)
end

@recipe function f(Xdom::DomainType{N, T}) where {N, T}
    opacity := 0.2
    color := :yellow
    legend := false

    dims = [1, 2]
    grid = get_grid(Xdom)
    dict = Dict{NTuple{2, Int}, Any}()
    for pos in enum_pos(Xdom)
        if !haskey(dict, pos[dims])
            dict[pos[dims]] = true
            @series begin
                return grid, pos
            end
        end
    end
end
