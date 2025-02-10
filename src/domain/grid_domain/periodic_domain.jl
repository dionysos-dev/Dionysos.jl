"""
    PeriodicDomainList{N, T, S} <: GridDomainType{N, T}

A periodic extension of `DomainList`, where specified dimensions wrap around at given periods.

# Fields
- `domain::DomainList{N, T, S}` : The underlying domain structure.
- `periodic_dims::Vector{Bool}` : Boolean vector indicating which dimensions are periodic.
- `periods::SVector{N, T}` : The length of the period in each dimension.
- `start::SVector{N, T}` : The start of the periodic cycle in each dimension.

# Features
- Uses `DomainList` for element management.
- Implements periodic boundary conditions for wrapping elements.
"""
struct PeriodicDomainList{N, T, S <: Grid{N, T}} <: GridDomainType{N, T}
    underlying_domain::DomainList{N, T, S}
    periodic_dims::Vector{Bool}
    periods::SVector{N, T}   
    start::SVector{N, T}
end

"""
    PeriodicDomainList(grid::S, periodic_dims::Vector{Bool}, periods::SVector{N, T}, start::SVector{N, T}) where {N, T, S <: Grid{N, T}}

Creates a periodic domain list with the given periodicity settings.
"""
function PeriodicDomainList(grid::S, periodic_dims::Vector{Bool}, periods::SVector{N, T}, start::SVector{N, T}) where {N, T, S <: Grid{N, T}}
    underlying_domain = DomainList(grid)
    return PeriodicDomainList(underlying_domain, periodic_dims, periods, start)
end

# ----------------------------
# Periodic Wrapping Logic
# ----------------------------

get_periodic_dims(domain::PeriodicDomainList) = domain.periodic_dims
get_period_lengths(domain::PeriodicDomainList) = domain.periods
get_period_start(domain::PeriodicDomainList) = domain.start

"""
    has_same_periodicity(domain1::PeriodicDomainList, domain2::PeriodicDomainList) -> Bool

Returns `true` if both domains have the same periodic settings.
"""
function has_same_periodicity(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    return domain1.periodic_dims == domain2.periodic_dims &&
           domain1.periods == domain2.periods &&
           domain1.start == domain2.start
end

"""
    wrap_pos(domain::PeriodicDomainList, pos::NTuple{N, Int}) -> NTuple{N, Int}

Transforms a position into its wrapped equivalent within the periodic boundaries.
"""
function wrap_pos(domain::PeriodicDomainList, pos::NTuple{N, Int})
    return Tuple(
        domain.periodic_dims[d] ? mod(pos[d] - domain.start[d], domain.periods[d]) + domain.start[d] : pos[d]
        for d in 1:length(pos)
    )
end


# ----------------------------
#  Overriding DomainList Methods for Periodicity
# ----------------------------

get_grid(domain::PeriodicDomainList) = get_grid(domain.underlying_domain)
enum_pos(domain::PeriodicDomainList) = enum_pos(domain.underlying_domain)
get_ncells(domain::PeriodicDomainList) = get_ncells(domain.underlying_domain) 
get_somepos(domain::PeriodicDomainList) = get_somepos(domain.underlying_domain)
Base.isempty(domain::PeriodicDomainList) = isempty(domain.underlying_domain)
Base.in(pos, domain::PeriodicDomainList) = in(wrap_pos(domain, pos), domain.underlying_domain)

get_pos_by_coord(domain::PeriodicDomainList, coord) = wrap_pos(domain, get_pos_by_coord(domain.underlying_domain, coord)) 
get_coord_by_pos(domain::PeriodicDomainList, pos) = get_coord_by_pos(domain.underlying_domain, wrap_pos(domain, pos))

function Base.issubset(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    if !has_same_periodicity(domain1, domain2)
        error("Cannot perform issubset on domains with different periodic settings.")
    end
    return issubset(domain1.underlying_domain, domain2.underlying_domain)
end

add_pos!(domain::PeriodicDomainList, pos) = add_pos!(domain.underlying_domain, wrap_pos(domain, pos))

function Base.union!(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    if !has_same_periodicity(domain1, domain2)
        error("Cannot perform union! on domains with different periodic settings.")
    end
    return union!(domain1.underlying_domain, domain2.underlying_domain)
end

function Base.setdiff!(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    if !has_same_periodicity(domain1, domain2)
        error("Cannot perform setdiff! on domains with different periodic settings.")
    end
    return setdiff!(domain1.underlying_domain, domain2.underlying_domain)
end

Base.empty!(domain::PeriodicDomainList) = empty!(domain.underlying_domain)
remove_pos!(domain::PeriodicDomainList, pos) = remove_pos!(domain.underlying_domain, wrap_pos(domain, pos))
