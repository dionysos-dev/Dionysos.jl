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

Creates a periodic domain list using an **existing grid**.

Ensures that:
1. The grid origin (`orig`) aligns with periodicity (`start + h/2` for periodic dimensions).
2. The grid step `h` is a multiple of the period in periodic dimensions.

Throws an error if misalignment occurs.
"""
function PeriodicDomainList(
    grid::S,
    periodic_dims::Vector{Bool},
    periods::SVector{N, T},
    start::SVector{N, T},
) where {N, T, S <: Grid{N, T}}
    orig = get_origin(grid)
    h = get_h(grid)

    for d in 1:N
        if periodic_dims[d]
            expected_orig = start[d] + h[d] / 2.0  # Expected alignment for periodic dimension

            # Ensure the grid origin aligns correctly
            if !isapprox(orig[d], expected_orig; atol = 1e-9)
                error(
                    "Grid origin orig[$d] = $(orig[d]) must be at start[$d] + h[$d]/2 = $(expected_orig) in periodic dimension $d.",
                )
            end

            # Ensure the period is a multiple of the grid step size
            if !isapprox(mod(periods[d], h[d]), 0.0; atol = 1e-9)
                error(
                    "Grid step size h[$d] = $(h[d]) must be a multiple of the period[$d] = $(periods[d]) in periodic dimension $d.",
                )
            end
        end
    end

    # Create the domain list and return the periodic domain
    underlying_domain = DomainList(grid)
    return PeriodicDomainList(underlying_domain, periodic_dims, periods, start)
end

"""
    PeriodicDomainList(periodic_dims::Vector{Bool}, periods::SVector{N, T}, start::SVector{N, T}, h::SVector{N, T}) where {N, T}

Creates a periodic domain list by automatically defining a grid with:
1. The origin aligns correctly with periodicity (`start + h/2` for periodic dimensions).
2. The grid step `h` is a multiple of `periods` in periodic dimensions.

Ensures that the **origin aligns with `start` and `periods`** in periodic dimensions.
Throws an error if misalignment occurs.
"""
function PeriodicDomainList(
    periodic_dims::Vector{Bool},
    periods::SVector{N, T},
    start::SVector{N, T},
    h::SVector{N, T},
) where {N, T}
    # Compute origin based on periodicity
    orig = SVector{N, T}([periodic_dims[d] ? start[d] + h[d] / 2.0 : start[d] for d in 1:N])

    # Ensure grid alignment with start and periods
    for d in 1:N
        if periodic_dims[d]
            if !isapprox(mod(periods[d], h[d]), 0.0; atol = 1e-9)
                error(
                    "Grid step size h[$d] = $(h[d]) must be a multiple of the period[$d] = $(periods[d]) in periodic dimension $d.",
                )
            end
        end
    end

    # Create the grid and return the periodic domain list
    grid = GridFree(orig, h)
    return PeriodicDomainList(grid, periodic_dims, periods, start)
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

Wraps a grid position `pos` into its equivalent within periodic boundaries.

**Handles periodicity correctly based on how the grid is defined.**
"""
function wrap_pos(domain::PeriodicDomainList{N, T}, pos::NTuple{N, Int}) where {N, T}
    h = get_h(get_grid(domain))
    start = domain.start
    periods = domain.periods

    return Tuple(
        domain.periodic_dims[d] ?
        mod(
            pos[d] - round(Int, (start[d] + h[d] / 2) / h[d]),
            round(Int, periods[d] / h[d]),
        ) + round(Int, (start[d] + h[d] / 2) / h[d]) : pos[d] for d in 1:N
    )
end

"""
    wrap_coord(domain::PeriodicDomainList, coord::SVector{N, Float64}) -> SVector{N, Float64}

Transforms a coordinate into its wrapped equivalent within the periodic boundaries.
"""
function wrap_coord(domain::PeriodicDomainList, coord::SVector{N, Float64}) where {N}
    return SVector(
        domain.periodic_dims[d] ?
        mod(coord[d] - domain.start[d], domain.periods[d]) + domain.start[d] : coord[d] for
        d in 1:N
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
        pos in Iterators.product(_ranges(rectI)...) if pos ∈ domain
    ]
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
