"""
    PeriodicGridMapping{N, T, M, P} <: GridMapping{N, T}

A periodic wrapper over an underlying `GridMapping`.

Fields (mirrors PeriodicDomainList):
- `periodic_dims::SVector{P, Int}`
- `periods::SVector{P, T}`
- `start::SVector{P, T}`
- `underlying_mapping::M`
- `periodic_index_map::NTuple{N, Union{Nothing, Int}}`
"""
struct PeriodicGridMapping{N, T, M <: GridMapping{N, T}, P} <: GridMapping{N, T}
    periodic_dims::SVector{P, Int}
    periods::SVector{P, T}
    start::SVector{P, T}
    underlying_mapping::M
    periodic_index_map::NTuple{N, Union{Nothing, Int}}
end

function _make_periodic_index_map(periodic_dims::SVector{P, Int}, N::Int) where {P}
    return ntuple(d -> begin
        i = findfirst(isequal(d), periodic_dims)
        isnothing(i) ? nothing : i
    end, N)
end

function get_grid_in_periods(
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
    return GridFree(SVector{N}(orig), h)
end

function get_grid_in_periods(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    h::SVector{N, T},
) where {N, T, P}
    start = zeros(SVector{P, T})
    return get_grid_in_periods(periodic_dims, periods, start, h)
end

# ----------------------------
# Constructors 
# ----------------------------

# NOTE (v0.2): two constructors taking a raw `grid::Grid` or a bare `h::SVector` were removed. Both
# were broken (one referenced an out-of-scope `mapping`; the other used a type parameter `M` absent
# from its `where` clause) and unused — every caller builds a `PeriodicGridMapping` from a
# `GridMapping` via the constructor below. Add a build-from-grid path, with tests, if it is needed.

function PeriodicGridMapping(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    mapping::M,
) where {N, T, M <: GridMapping{N, T}, P}
    grid = get_grid(mapping)
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

    pmap = _make_periodic_index_map(periodic_dims, N)
    return PeriodicGridMapping{N, T, M, P}(periodic_dims, periods, start, mapping, pmap)
end

function PeriodicGridMapping(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    mapping::M,
) where {N, T, M <: GridMapping{N, T}, P}
    start = zeros(SVector{P, T})
    return PeriodicGridMapping(periodic_dims, periods, start, mapping)
end

# Non-periodic wrapper (convenience), same as PeriodicDomainList(grid::Grid).
function PeriodicGridMapping(mapping::M) where {N, T, M <: GridMapping{N, T}}
    return PeriodicGridMapping(
        SVector{0, Int}(),
        SVector{0, T}(),
        zeros(SVector{0, T}),
        mapping,
    )
end

# ----------------------------
# Periodic metadata API
# ----------------------------

get_periodic_dims(m::PeriodicGridMapping) = m.periodic_dims
get_periods(m::PeriodicGridMapping) = m.periods
get_periodic_starts(m::PeriodicGridMapping) = m.start

is_periodic(m::PeriodicGridMapping, d::Int) = m.periodic_index_map[d] !== nothing
is_periodic(m::PeriodicGridMapping) =
    any(d -> is_periodic(m, d), 1:length(m.periodic_index_map))

function has_same_periodicity(m1::PeriodicGridMapping, m2::PeriodicGridMapping)
    return get_periodic_dims(m1) == get_periodic_dims(m2) &&
           get_periods(m1) == get_periods(m2) &&
           get_periodic_starts(m1) == get_periodic_starts(m2)
end

# ----------------------------
# Wrapping logic (pos + coord)
# ----------------------------

function wrap_pos(m::PeriodicGridMapping{N, T}, pos::NTuple{N, Int}) where {N, T}
    !is_periodic(m) && return pos

    h = get_h(get_grid(m))
    pmap = m.periodic_index_map
    periods = m.periods

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

function wrap_coord(m::PeriodicGridMapping{N, T}, x::SVector{N, T}) where {N, T}
    !is_periodic(m) && return x

    pmap = m.periodic_index_map
    return SVector{N, T}(
        ntuple(
            d -> begin
                i = pmap[d]
                i === nothing ? x[d] : UT.wrap_value(x[d], m.start[i], m.periods[i])
            end,
            N,
        ),
    )
end

# ----------------------------
# Override GridMapping methods by delegation + wrapping
# ----------------------------

get_grid(m::PeriodicGridMapping) = get_grid(m.underlying_mapping)

get_n_state(m::PeriodicGridMapping) = get_n_state(m.underlying_mapping)
enum_states(m::PeriodicGridMapping) = enum_states(m.underlying_mapping)

get_pos_by_state(m::PeriodicGridMapping{N}, q::Int) where {N} =
    get_pos_by_state(m.underlying_mapping, q)

# validity in pos-space: wrap first, then check underlying
is_valid_pos(m::PeriodicGridMapping{N}, pos::NTuple{N, Int}) where {N} =
    is_valid_pos(m.underlying_mapping, wrap_pos(m, pos))

# pos -> state: wrap first, then delegate
get_state_by_pos(m::PeriodicGridMapping{N}, pos::NTuple{N, Int}) where {N} =
    get_state_by_pos(m.underlying_mapping, wrap_pos(m, pos))

# coord -> state: wrap coordinate first (like your domain wrap_coord), then do normal path
get_state_by_coord(m::PeriodicGridMapping, x) =
    get_state_by_coord(m.underlying_mapping, wrap_coord(m, x))

# coord by state can remain the base (it uses underlying pos by state)
get_coord_by_state(m::PeriodicGridMapping, q::Int) =
    get_coord_by_pos(get_grid(m), get_pos_by_state(m, q))

# Only valid when underlying mappig is Explicit
add_pos!(m::PeriodicGridMapping, pos) = add_pos!(m.underlying_mapping, wrap_pos(m, pos))
function add_set!(m::PeriodicGridMapping, set, incl_mode::INCL_MODE)
    grid = get_grid(m)
    for pos in get_pos_from_set(grid, set, incl_mode)
        add_pos!(m, pos)
    end
    return m
end
