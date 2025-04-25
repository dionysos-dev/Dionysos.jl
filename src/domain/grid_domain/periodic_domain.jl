"""
    PeriodicDomainList{N, T, S, P} <: GridDomainType{N, T}

A periodic extension of `DomainList`, where specified dimensions wrap around with specified periods.

# Fields
- `periodic_dims::SVector{P, Int}`: Indices of the periodic dimensions.
- `periods::SVector{P, T}`: Period length in each periodic dimension.
- `start::SVector{P, T}`: Start value for each periodic dimension.
- `underlying_domain::DomainList{N, T, S}`: The wrapped domain over all `N` dimensions.
"""
struct PeriodicDomainList{N, T, S <: Grid{N, T}, P} <: GridDomainType{N, T}
    periodic_dims::SVector{P, Int}
    periods::SVector{P, T}
    start::SVector{P, T}
    underlying_domain::DomainList{N, T, S}
    periodic_index_map::NTuple{N, Union{Nothing, Int}}
end

"""
    _make_periodic_index_map(periodic_dims::SVector{P, Int}, N::Int)

Returns an `NTuple{N, Union{Nothing, Int}}` where each entry is either:
- `nothing` if dimension `d` is not periodic
- `i` such that `periodic_dims[i] == d`
"""
function _make_periodic_index_map(periodic_dims::SVector{P, Int}, N::Int) where {P}
    return ntuple(d -> begin
        i = findfirst(isequal(d), periodic_dims)
        isnothing(i) ? nothing : i
    end, N)
end

"""
    PeriodicDomainList(
        periodic_dims::SVector{P, Int},
        periods::SVector{P, T},
        start::SVector{P, T},
        grid::S
    ) where {N, T, S <: Grid{N, T}, P}

Constructs a `PeriodicDomainList` from an existing grid.

# Requirements
For each periodic dimension `d = periodic_dims[i]`:
- The grid origin must satisfy `origin[d] == start[i] + h[d] / 2`
- The grid step `h[d]` must divide the period `periods[i]`

Throws an error if constraints are not met.
"""
function PeriodicDomainList(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    grid::S,
) where {N, T, S <: Grid{N, T}, P}
    orig = get_origin(grid)
    h = get_h(grid)

    for i in 1:P
        d = periodic_dims[i]
        expected_orig = start[i] + h[d] / 2.0

        if !isapprox(orig[d], expected_orig; atol = 1e-9)
            error(
                "Grid origin orig[$d] = $(orig[d]) must equal start[$i] + h[$d]/2 = $(expected_orig).",
            )
        end

        q = periods[i] / h[d]
        if !isapprox(q, round(q); atol = 1e-9)
            error("Grid step h[$d] = $(h[d]) must divide period[$i] = $(periods[i]).")
        end
    end
    domain = DomainList(grid)
    map = _make_periodic_index_map(periodic_dims, N)
    return PeriodicDomainList{N, T, S, P}(periodic_dims, periods, start, domain, map)
end

function PeriodicDomainList(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    grid::S,
) where {N, T, S <: Grid{N, T}, P}
    start = zeros(SVector{P, T})
    return PeriodicDomainList(periodic_dims, periods, start, grid)
end

"""
    PeriodicDomainList(
        periodic_dims::SVector{P, Int},
        periods::SVector{P, T},
        start::SVector{P, T},
        h::SVector{N, T}
    ) where {N, T, P}

Constructs a periodic domain list by generating a grid whose origin aligns with the periodic structure:
- In each periodic dimension `d = periodic_dims[i]`, the origin is set to `start[i] + h[d]/2`
- All other dimensions default to origin = 0.0

Throws an error if the period is not divisible by the grid step.
"""
function PeriodicDomainList(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    h::SVector{N, T},
) where {N, T, P}
    orig = zeros(N)
    for i in 1:P
        d = periodic_dims[i]
        orig[d] = start[i] + h[d] / 2.0
    end

    for i in 1:P
        d = periodic_dims[i]
        q = periods[i] / h[d]
        if !isapprox(q, round(q); atol = 1e-9)
            error("Grid step h[$d] = $(h[d]) must divide period[$i] = $(periods[i]).")
        end
    end

    grid = GridFree(SVector{N}(orig), h)
    return PeriodicDomainList(periodic_dims, periods, start, grid)
end

"""
    PeriodicDomainList(
        periodic_dims::SVector{P, Int},
        periods::SVector{P, T},
        h::SVector{N, T}
    ) where {N, T, P}

Constructs a periodic domain list with `start = 0` in all periodic dimensions.
"""
function PeriodicDomainList(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    h::SVector{N, T},
) where {N, T, P}
    start = zeros(SVector{P, T})
    return PeriodicDomainList(periodic_dims, periods, start, h)
end

"""
    PeriodicDomainList(h::SVector{N, T}) where {N, T}

Creates a non-periodic `PeriodicDomainList` by wrapping a grid with uniform spacing `h`.

This is equivalent to a regular `DomainList` with periodicity disabled.
"""
function PeriodicDomainList(h::SVector{N, T}) where {N, T}
    return PeriodicDomainList(SVector{0, Int}(), SVector{0, T}(), h)
end

"""
    PeriodicDomainList(grid::Grid) -> PeriodicDomainList

Wraps a given `Grid` in a non-periodic `PeriodicDomainList`.

This is equivalent to a regular domain list with no periodic behavior applied.
"""
function PeriodicDomainList(grid::S) where {N, T, S <: Grid{N, T}}
    return PeriodicDomainList(SVector{0, Int}(), SVector{0, T}(), SVector{0, T}(), grid)
end

# ----------------------------
# Periodic Wrapping Logic
# ----------------------------

get_periodic_dims(domain::PeriodicDomainList) = domain.periodic_dims
get_periods(domain::PeriodicDomainList) = domain.periods
get_periodic_starts(domain::PeriodicDomainList) = domain.start

"""
    is_periodic(domain::PeriodicDomainList, d::Int) -> Bool

Returns `true` if dimension `d` is periodic in the given domain.
"""
is_periodic(domain::PeriodicDomainList, d::Int) = domain.periodic_index_map[d] !== nothing

"""
    is_periodic(domain::PeriodicDomainList) -> Bool

Returns `true` if **any** dimension in the domain is periodic.
Uses `is_periodic(domain, d)` internally.
"""
function is_periodic(domain::PeriodicDomainList)
    N = length(domain.periodic_index_map)
    return any(d -> is_periodic(domain, d), 1:N)
end

"""
    has_same_periodicity(domain1::PeriodicDomainList, domain2::PeriodicDomainList) -> Bool

Returns `true` if both domains have the same periodic settings.
"""
function has_same_periodicity(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    return get_periodic_dims(domain1) == get_periodic_dims(domain2) &&
           get_periods(domain1) == get_periods(domain2) &&
           get_periods(domain1) == get_periods(domain2)
end

"""
    wrap_pos(domain::PeriodicDomainList{N, T}, pos::NTuple{N, Int}) -> NTuple{N, Int}

Wraps a grid position `pos` into its equivalent within the periodic boundaries defined by `domain`.

Only the periodic dimensions (as defined in `domain.periodic_dims`) are wrapped.  
Non-periodic dimensions are returned as-is.

# Returns
- An `NTuple{N, Int}` where periodic dimensions are wrapped modulo the number of cells in the period.
"""
function wrap_pos(domain::PeriodicDomainList{N, T}, pos::NTuple{N, Int}) where {N, T}
    if !is_periodic(domain)
        return pos
    end

    h = get_h(get_grid(domain))
    start = get_periodic_starts(domain)
    periods = get_periods(domain)
    pmap = domain.periodic_index_map

    return ntuple(d -> begin
        i = pmap[d]
        if i === nothing
            pos[d]
        else
            span = round(Int, periods[i] / h[d])
            mod(pos[d], span)
        end
    end, N)
end

"""
    wrap_coord(domain::PeriodicDomainList{N, T}, coord::SVector{N, T}) -> SVector{N, T}

Wraps a real-valued coordinate `coord` into its periodic equivalent, based on the periodic
dimensions defined in the `domain`.

Non-periodic dimensions are returned unchanged.
"""
function wrap_coord(domain::PeriodicDomainList{N, T}, coord::SVector{N, T}) where {N, T}
    if !is_periodic(domain)
        return coord
    end

    pmap = domain.periodic_index_map
    return SVector{N, T}(ntuple(d -> begin
        i = pmap[d]
        if i === nothing
            coord[d]
        else
            s = domain.start[i]
            p = domain.periods[i]
            mod(coord[d] - s, p) + s
        end
    end, N))
end

# ----------------------------
#  Overriding GridDomainType Methods for Periodicity
# ----------------------------

get_grid(domain::PeriodicDomainList) = get_grid(domain.underlying_domain)
enum_pos(domain::PeriodicDomainList) = enum_pos(domain.underlying_domain)
get_ncells(domain::PeriodicDomainList) = get_ncells(domain.underlying_domain)
get_somepos(domain::PeriodicDomainList) = get_somepos(domain.underlying_domain)
Base.isempty(domain::PeriodicDomainList) = isempty(domain.underlying_domain)
Base.in(pos, domain::PeriodicDomainList) =
    in(wrap_pos(domain, pos), domain.underlying_domain)

get_pos_by_coord(domain::PeriodicDomainList, coord) =
    wrap_pos(domain, get_pos_by_coord(domain.underlying_domain, coord))
get_coord_by_pos(domain::PeriodicDomainList, pos) =
    get_coord_by_pos(domain.underlying_domain, wrap_pos(domain, pos))
function get_subset_pos(
    domain::PeriodicDomainList,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(get_grid(domain), rect, incl_mode)
    return [
        wrap_pos(domain, pos) for
        pos in Iterators.product(_ranges(rectI)...) if in(pos, domain)
    ]
end

function get_subset_pos_in_grid(
    domain::PeriodicDomainList,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(get_grid(domain), rect, incl_mode)
    raw_iter = Iterators.product(_ranges(rectI)...)
    return (wrap_pos(domain, pos) for pos in raw_iter)
end

function Base.issubset(domain1::PeriodicDomainList, domain2::PeriodicDomainList)
    if !has_same_periodicity(domain1, domain2)
        error("Cannot perform issubset on domains with different periodic settings.")
    end
    return issubset(domain1.underlying_domain, domain2.underlying_domain)
end

add_pos!(domain::PeriodicDomainList, pos) =
    add_pos!(domain.underlying_domain, wrap_pos(domain, pos))

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
remove_pos!(domain::PeriodicDomainList, pos) =
    remove_pos!(domain.underlying_domain, wrap_pos(domain, pos))
