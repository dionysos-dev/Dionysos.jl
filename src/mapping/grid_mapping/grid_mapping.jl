"""
GridMapping{N,T} <: AbstractMapping{N,T}:
- adds grid, and pos <-> state conversions
"""
abstract type GridMapping{N,T} <: AbstractMapping{N,T} end

get_grid(m::GridMapping) = error("implement get_grid")
get_state_by_pos(m::GridMapping{N}, pos::NTuple{N,Int}) where {N} =
    error("implement get_state_by_pos")
get_pos_by_state(m::GridMapping{N}, q::Int) where {N} =
    error("implement get_pos_by_state")
is_valid_pos(m::GridMapping{N}, pos::NTuple{N,Int}) where {N} = true


enum_pos(m::GridMapping{N,T}) where {N,T} = (get_pos_by_state(m, q) for q in enum_states(m))
get_elem_by_pos(m::GridMapping{N,T}, pos::NTuple{N,Int}) where {N,T} = get_elem_by_pos(get_grid(m), pos)

get_elem_by_state(m::GridMapping{N,T}, q::Int) where {N,T} = get_elem_by_pos(m, get_pos_by_state(m, q))
get_elem_by_coord(m::GridMapping{N,T}, x) where {N,T} = get_elem_by_pos(m, get_pos_by_coord(get_grid(m), x))
get_state_by_coord(m::GridMapping, x) = get_state_by_pos(m, get_pos_by_coord(get_grid(m), x))
get_coord_by_state(m::GridMapping{N,T}, q::Int) where {N,T} = get_coord_by_pos(get_grid(m), get_pos_by_state(m, q))

convert_to_list_mapping(m::GridMapping) = ListMapping(collect(enum_coords(m)))


# returns (states, allin) where allin=false if any covered grid-pos is invalid
function get_states_from_set_strict(
    m::AbstractMapping{N,T},
    x::AbstractVector,
    ::INCL_MODE,
) where {N,T}
    q = get_state_by_coord(m, x)
    return q===nothing ? (nothing, false) : (Int[q], true)
end

function get_states_from_set_strict(
    m::GridMapping{N},
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
) where {N}
    grid = get_grid(m)
    rectI = get_pos_lims(grid, rect, incl_mode)

    qs = Int[]
    allin = true

    for pos in Iterators.product(_ranges(rectI)...)
        p = pos::NTuple{N,Int}
        if is_valid_pos(m, p)
            push!(qs, get_state_by_pos(m, p))
        else
            allin = false
        end
    end
    return qs, allin
end

function get_states_from_set_strict(
    m::GridMapping{N},
    subsets::UT.LazyUnion,
    incl_mode::INCL_MODE,
) where {N}
    acc = Int[]
    allin = true
    for subset in subsets.sets
        qs, ok = get_states_from_set_strict(m, subset, incl_mode)
        Base.append!(acc, qs)
        allin &= ok
    end
    unique!(acc)
    return acc, allin
end

function get_states_from_set_strict(
    m::GridMapping{N},
    set::UT.LazySetMinus,
    incl_mode::INCL_MODE,
) where {N}
    states_A, okA = get_states_from_set_strict(m, set.A, incl_mode)
    states_B, okB = get_states_from_set_strict(m, set.B, _invInclMode(incl_mode))
    return setdiff(states_A, states_B), (okA && okB)
end

function get_states_from_set(
    m::AbstractMapping{N,T},
    x::AbstractVector,
    incl_mode::INCL_MODE,
) where {N,T}
    qs, _ = get_states_from_set_strict(m, x, incl_mode)
    return qs
end

function get_states_from_set(
    m::GridMapping{N},
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
) where {N}
    qs, _ = get_states_from_set_strict(m, rect, incl_mode)
    return qs
end

function get_states_from_set(
    m::GridMapping{N},
    subsets::UT.LazyUnion,
    incl_mode::INCL_MODE,
) where {N}
    qs, _ = get_states_from_set_strict(m, subsets, incl_mode)
    return qs
end

function get_states_from_set(
    m::GridMapping{N},
    set::UT.LazySetMinus,
    incl_mode::INCL_MODE,
) where {N}
    qs, _ = get_states_from_set_strict(m, set, incl_mode)
    return qs
end

# ----------------------------
# Plotting
# ----------------------------

struct IntRect2
    lb::NTuple{2,Int}
    ub::NTuple{2,Int}
end

function merge_rectangles_2d(pos2d::AbstractVector{NTuple{2,Int}})
    isempty(pos2d) && return IntRect2[]

    rows = Dict{Int, Vector{Int}}()
    for (x,y) in pos2d
        push!(get!(rows, y, Int[]), x)
    end
    for xs in values(rows); sort!(xs) end

    intervals = Tuple{Int,Int,Int}[]  # (y, x1, x2)
    ys = sort!(collect(keys(rows)))
    for y in ys
        xs = rows[y]
        i = 1
        while i <= length(xs)
            x1 = xs[i]; x2 = x1
            while i < length(xs) && xs[i+1] == x2 + 1
                i += 1; x2 = xs[i]
            end
            push!(intervals, (y, x1, x2))
            i += 1
        end
    end

    sort!(intervals)
    rects = IntRect2[]
    i = 1
    while i <= length(intervals)
        y, x1, x2 = intervals[i]
        y1 = y; y2 = y
        i += 1
        while i <= length(intervals)
            y_next, x1_next, x2_next = intervals[i]
            if x1_next == x1 && x2_next == x2 && y_next == y2 + 1
                y2 = y_next
                i += 1
            else
                break
            end
        end
        push!(rects, IntRect2((x1,y1), (x2,y2)))
    end
    return rects
end

function intrect2_to_real_rect(grid, r::IntRect2, d1::Int, d2::Int)
    orig = get_origin(grid)
    h    = get_h(grid)
    T    = eltype(orig)

    lb = SVector{2,T}(
        orig[d1] + r.lb[1]*h[d1] - h[d1]/2,
        orig[d2] + r.lb[2]*h[d2] - h[d2]/2,
    )
    ub = SVector{2,T}(
        orig[d1] + r.ub[1]*h[d1] + h[d1]/2,
        orig[d2] + r.ub[2]*h[d2] + h[d2]/2,
    )
    return UT.HyperRectangle(lb, ub)
end

"""function positions_from_states_unique_on_dims(
    m::GridMapping{N},
    states;
    dims = [1, 2],
) where {N}
    d1, d2 = dims[1], dims[2]
    out  = NTuple{N,Int}[]
    seen = Set{NTuple{2,Int}}()

    for q in states
        p = get_pos_by_state(m, q) 
        key = (p[d1], p[d2]) 
        if !(key in seen)
            push!(seen, key)
            push!(out, p)
        end
    end
    return out
end"""

"""
Project states onto (d1,d2). For each 2D cell keep:
- best value (min by default)
- one representative full pos (NTuple{N,Int}) to plot the cell (slow branch),
  or just to build pos2d list (fast branch).
"""
function project_states_on_dims(
    m::GridMapping{N},
    states;
    dims = [1,2],
    value_function::Union{Nothing,Function} = nothing,
    reducer::Function = min,
) where {N}
    d1, d2 = dims[1], dims[2]

    if value_function === nothing
        # just dedup by (d1,d2), keep any representative full pos
        seen = Set{NTuple{2,Int}}()
        posN = NTuple{N,Int}[]
        for q in states
            p = get_pos_by_state(m, q)
            key = (p[d1], p[d2])
            if !(key in seen)
                push!(seen, key)
                push!(posN, p)
            end
        end
        return nothing, posN
    else
        # map: key -> (best_value, representative_full_pos)
        proj = Dict{NTuple{2,Int}, Tuple{Float64, NTuple{N,Int}}}()
        for q in states
            p = get_pos_by_state(m, q)
            key = (p[d1], p[d2])
            v = Float64(value_function(q))
            if haskey(proj, key)
                vbest, _ = proj[key]
                newbest = reducer(vbest, v)
                if newbest != vbest
                    proj[key] = (newbest, p)
                end
            else
                proj[key] = (v, p)
            end
        end
        posN = [vp[2] for vp in values(proj)]
        return proj, posN
    end
end

# ----------------------------
# Concrete grid mappings
# ----------------------------

"""
ExplicitGridMapping:
explicit enumeration of positions -> ids (dict) and ids -> positions (vector)
"""
struct ExplicitGridMapping{N,T,G} <: GridMapping{N,T}
    grid::G
    pos2id::Dict{NTuple{N,Int},Int}
    id2pos::Vector{NTuple{N,Int}}
end

function ExplicitGridMapping(grid::G) where {N,T,G}
    return ExplicitGridMapping{N,T,G}(
        grid,
        Dict{NTuple{N,Int},Int}(),
        NTuple{N,Int}[],
    )
end

function ExplicitGridMapping{N,T}(grid::G) where {N,T,G}
    return ExplicitGridMapping{N,T,G}(
        grid,
        Dict{NTuple{N,Int},Int}(),
        NTuple{N,Int}[],
    )
end

function ExplicitGridMapping(grid::Grid{N,T}) where {N,T}
    return ExplicitGridMapping{N,T}(grid)
end

function ExplicitGridMapping{N,T}(grid::G, positions) where {N,T,G}
    m = ExplicitGridMapping{N,T}(grid)
    for pos in positions
        add_pos!(m, pos)   # dedup via pos2id
    end
    return m
end

function ExplicitGridMapping{N,T}(grid::G, set, incl_mode::INCL_MODE) where {N,T,G}
    m = ExplicitGridMapping{N,T}(grid)
    add_set!(m, set, incl_mode)
    return m
end


get_grid(m::ExplicitGridMapping) = m.grid
get_n_state(m::ExplicitGridMapping) = length(m.id2pos)
enum_states(m::ExplicitGridMapping) = 1:get_n_state(m)

is_valid_pos(m::ExplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} = haskey(m.pos2id, pos)
get_state_by_pos(m::ExplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} = is_valid_pos(m, pos) ? m.pos2id[pos] : nothing
get_pos_by_state(m::ExplicitGridMapping{N}, q::Int) where {N} = m.id2pos[q]

function add_pos!(m::ExplicitGridMapping{N,T}, pos::NTuple{N,Int}) where {N,T}
    if haskey(m.pos2id, pos)
        return m.pos2id[pos]
    end
    q = length(m.id2pos) + 1
    push!(m.id2pos, pos)
    m.pos2id[pos] = q
    return q
end

function add_set!(m::ExplicitGridMapping{N,T}, set, incl_mode::INCL_MODE) where {N,T}
    grid = get_grid(m)
    for pos in get_pos_from_set(grid, set, incl_mode)
        add_pos!(m, pos)   # dedup happens here
    end
    return m
end


"""
ImplicitGridMapping:
- universe is a rectangular index box (min_pos..max_pos) possibly periodic
- id computed by row-major linearization (no dict)
- pos recovered by inverse mapping
"""
struct ImplicitGridMapping{N,T,G} <: GridMapping{N,T}
    grid::G
    min_pos::NTuple{N,Int}
    max_pos::NTuple{N,Int}
    n_per_dim::NTuple{N,Int}
    strides::NTuple{N,Int}
end

function ImplicitGridMapping(
    grid::G,
    min_pos::NTuple{N,Int},
    max_pos::NTuple{N,Int},
) where {N,G}
    n_per_dim = ntuple(d -> max_pos[d] - min_pos[d] + 1, N)
    @assert all(n_per_dim .>= 1) "Invalid implicit box: some max_pos < min_pos"

    # row-major strides: stride[1]=1, stride[2]=n1, stride[3]=n1*n2, ...
    strides = ntuple(d -> begin
        s = 1
        for k in 1:(d-1)
            s *= n_per_dim[k]
        end
        s
    end, N)

    T = eltype(get_origin(grid))
    return ImplicitGridMapping{N,T,G}(grid, min_pos, max_pos, n_per_dim, strides)
end

function ImplicitGridMapping(
    grid::G,
    rect::UT.HyperRectangle{N,T};
    incl_mode=OUTER
) where {N,T,G}
    rectI = get_pos_lims(grid, rect, incl_mode)
    min_pos = rectI.lb isa NTuple{N,Int} ? rectI.lb : NTuple{N,Int}(rectI.lb)
    max_pos = rectI.ub isa NTuple{N,Int} ? rectI.ub : NTuple{N,Int}(rectI.ub)
    return ImplicitGridMapping(grid, min_pos, max_pos)
end

get_grid(m::ImplicitGridMapping) = m.grid
get_n_state(m::ImplicitGridMapping) = prod(m.n_per_dim)
enum_states(m::ImplicitGridMapping) = 1:get_n_state(m)

is_valid_pos(m::ImplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} =
    all(d -> (m.min_pos[d] <= pos[d] <= m.max_pos[d]), 1:N)

function get_state_by_pos(m::ImplicitGridMapping{N}, pos::NTuple{N,Int}) where {N}
    is_valid_pos(m, pos) || throw(DomainError(pos, "pos out of implicit mapping"))
    q = 1
    @inbounds for d in 1:N
        q += (pos[d] - m.min_pos[d]) * m.strides[d]
    end
    return q
end

function get_pos_by_state(m::ImplicitGridMapping{N}, q::Int) where {N}
    (1 <= q <= get_n_state(m)) || throw(DomainError(q, "state out of implicit universe"))
    z = q - 1
    return ntuple(d -> begin
        span  = m.n_per_dim[d]
        digit = (z ÷ m.strides[d]) % span
        m.min_pos[d] + digit
    end, N)
end



# ----------------------------
# Periodic wrapper (decorator)
# ----------------------------

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
struct PeriodicGridMapping{N, T, M <: GridMapping{N,T}, P} <: GridMapping{N,T}
    periodic_dims::SVector{P, Int}
    periods::SVector{P, T}
    start::SVector{P, T}
    underlying_mapping::M
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
Constructs a grid whose origin aligns with the periodic structure:
- In each periodic dimension `d = periodic_dims[i]`, the origin is set to `start[i] + h[d]/2`
- All other dimensions default to origin = 0.0

Throws an error if the period is not divisible by the grid step.
"""
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

"""
    PeriodicGridMapping(periodic_dims, periods, start, mapping)

Checks:
- origin[d] == start[i] + h[d]/2
- h[d] divides periods[i]
"""
function PeriodicGridMapping(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    grid::Grid{N,T},
) where {N, T, M <: GridMapping{N,T}, P}

    orig = get_origin(grid)
    h = get_h(grid)

    for i in 1:P
        d = periodic_dims[i]
        expected_orig = start[i] + h[d] / 2.0
        if !isapprox(orig[d], expected_orig; atol = 1e-9)
            error("Grid origin orig[$d] = $(orig[d]) must equal start[$i] + h[$d]/2 = $(expected_orig).")
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
    start::SVector{P, T},
    h::SVector{N, T},
) where {N, T, P}

    grid = get_grid_in_periods(periodic_dims, periods, start, h)
    return PeriodicGridMapping{N, T, M, P}(periodic_dims, periods, start, grid)
end


"""
    PeriodicGridMapping(periodic_dims, periods, start, mapping)

Checks:
- origin[d] == start[i] + h[d]/2
- h[d] divides periods[i]
"""
function PeriodicGridMapping(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
    mapping::M,
) where {N, T, M <: GridMapping{N,T}, P}

    grid = get_grid(mapping)
    orig = get_origin(grid)
    h = get_h(grid)

    for i in 1:P
        d = periodic_dims[i]
        expected_orig = start[i] + h[d] / 2.0
        if !isapprox(orig[d], expected_orig; atol = 1e-9)
            error("Grid origin orig[$d] = $(orig[d]) must equal start[$i] + h[$d]/2 = $(expected_orig).")
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
) where {N, T, M <: GridMapping{N,T}, P}
    start = zeros(SVector{P, T})
    return PeriodicGridMapping(periodic_dims, periods, start, mapping)
end

"""
Non-periodic wrapper (convenience), same as PeriodicDomainList(grid::Grid).
"""
function PeriodicGridMapping(mapping::M) where {N, T, M <: GridMapping{N,T}}
    return PeriodicGridMapping(SVector{0, Int}(), SVector{0, T}(), zeros(SVector{0, T}), mapping)
end

# ----------------------------
# Periodic metadata API
# ----------------------------

get_periodic_dims(m::PeriodicGridMapping) = m.periodic_dims
get_periods(m::PeriodicGridMapping) = m.periods
get_periodic_starts(m::PeriodicGridMapping) = m.start

is_periodic(m::PeriodicGridMapping, d::Int) = m.periodic_index_map[d] !== nothing
is_periodic(m::PeriodicGridMapping) = any(d -> is_periodic(m, d), 1:length(m.periodic_index_map))

function has_same_periodicity(m1::PeriodicGridMapping, m2::PeriodicGridMapping)
    return get_periodic_dims(m1) == get_periodic_dims(m2) &&
           get_periods(m1) == get_periods(m2) &&
           get_periodic_starts(m1) == get_periodic_starts(m2)
end

# ----------------------------
# Wrapping logic (pos + coord)
# ----------------------------

"""
    wrap_pos(m, pos)

Same convention as your PeriodicDomainList:
span = round(Int, periods[i]/h[d])
wrapped = mod(pos[d], span)

If your grid pos are 1-based, replace with:
    1 + mod(pos[d]-1, span)
"""
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
    return SVector{N, T}(ntuple(d -> begin
        i = pmap[d]
        if i === nothing
            x[d]
        else
            s = m.start[i]
            p = m.periods[i]
            mod(x[d] - s, p) + s
        end
    end, N))
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
is_valid_pos(m::PeriodicGridMapping{N}, pos::NTuple{N,Int}) where {N} =
    is_valid_pos(m.underlying_mapping, wrap_pos(m, pos))

# pos -> state: wrap first, then delegate
get_state_by_pos(m::PeriodicGridMapping{N}, pos::NTuple{N,Int}) where {N} =
    get_state_by_pos(m.underlying_mapping, wrap_pos(m, pos))

# coord -> state: wrap coordinate first (like your domain wrap_coord), then do normal path
get_state_by_coord(m::PeriodicGridMapping, x) = get_state_by_coord(m.underlying_mapping, wrap_coord(m, x))

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
