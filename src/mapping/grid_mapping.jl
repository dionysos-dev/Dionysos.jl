module Mappings

# ----------------------------
# Mapping API
# ----------------------------

export AbstractMapping, GridMapping,
       get_n_state, enum_states,
       get_state_by_coord, get_coord_by_state,
       get_grid,
       get_state_by_pos, get_pos_by_state,
       is_valid_state, is_valid_pos

"""
AbstractMapping{N,T}:
- defines a universe of labels (1..get_n_state)
- provides coord <-> state conversions
"""
abstract type AbstractMapping{N,T} end

get_n_state(m::AbstractMapping) = error("implement get_n_state(::AbstractMapping)")
enum_states(m::AbstractMapping) = 1:get_n_state(m)

get_state_by_coord(m::AbstractMapping, x) = error("implement get_state_by_coord")
get_coord_by_state(m::AbstractMapping, q::Int) = error("implement get_coord_by_state")

get_states_from_set(m::AbstractMapping, set, incl_mode) =
    error("implement get_states_from_set")

"Optional validity checks (default true for universe ids)."
is_valid_state(m::AbstractMapping, q::Int) = (1 <= q <= get_n_state(m))

"""
GridMapping{N,T} <: AbstractMapping{N,T}:
- adds grid, and pos <-> state conversions
"""
abstract type GridMapping{N,T} <: AbstractMapping{N,T} end

get_grid(m::GridMapping) = error("implement get_grid(::GridMapping)")

get_state_by_pos(m::GridMapping{N}, pos::NTuple{N,Int}) where {N} =
    error("implement get_state_by_pos(::GridMapping, pos)")

get_pos_by_state(m::GridMapping{N}, q::Int) where {N} =
    error("implement get_pos_by_state(::GridMapping, q)")

is_valid_pos(m::GridMapping{N}, pos::NTuple{N,Int}) where {N} = true

# Derived coord/state conversions for any GridMapping
get_state_by_coord(m::GridMapping, x) =
    get_state_by_pos(m, Dionysos.Domain.get_pos_by_coord(get_grid(m), x))

get_coord_by_state(m::GridMapping, q::Int) =
    Dionysos.Domain.get_coord_by_pos(get_grid(m), get_pos_by_state(m, q))

# ---------- Periodicity API (optional but recommended) ----------

"Return periodic dims flags (NTuple{N,Bool}). Default: no periodicity."
periodic_dims(m::GridMapping{N}) where {N} = ntuple(_ -> false, N)

"Wrap a pos into the mapping’s periodic box. Default: identity."
wrap_pos(m::GridMapping{N}, pos::NTuple{N,Int}) where {N} = pos

"Wrap a coordinate into the mapping’s periodic coordinate range. Default: identity."
wrap_coord(m::GridMapping, x) = x

function get_states_from_set(
    m::GridMapping{N},
    rect::UT.HyperRectangle,
    incl_mode::DO.INCL_MODE,
) where {N}
    grid = get_grid(m)
    rectI = DO.get_pos_lims(grid, rect, incl_mode)
    qs = Int[]
    # iterate all integer positions in that index-box
    for pos in Iterators.product(DO._ranges(rectI)...)
        p = pos::NTuple{N,Int}
        # Filter by mapping validity:
        # - explicit mapping: only enumerated positions
        # - implicit mapping: bounds check (with periodic wrap handled inside get_state_by_pos)
        if is_valid_pos(m, p)
            # For implicit periodic mappings: get_state_by_pos wraps.
            # For explicit: will KeyError if you don't filter; we filtered.
            push!(qs, get_state_by_pos(m, p))
        end
    end
    return qs
end

function get_states_from_set(
    m::GridMapping{N},
    subsets::UT.LazyUnionSetArray,
    incl_mode::DO.INCL_MODE,
) where {N}
    acc = Int[]
    for subset in subsets.sets
        append!(acc, get_states_from_set(m, subset, incl_mode))
    end
    unique!(acc)
    return acc
end

function get_states_from_set(
    m::GridMapping{N},
    set::UT.LazySetMinus,
    incl_mode::DO.INCL_MODE,
) where {N}
    # A uses incl_mode
    states_A = get_states_from_set(m, set.A, incl_mode)
    # B uses opposite incl mode (same convention you used before)
    states_B = get_states_from_set(m, set.B, DO._invInclMode(incl_mode))
    return setdiff(states_A, states_B)
end

# ----------------------------
# Concrete grid mappings
# ----------------------------

export ExplicitGridMapping, ImplicitGridMapping

"""
ExplicitGridMapping:
- explicit enumeration of positions -> ids (dict) and ids -> positions (vector)
- universe = exactly the enumerated positions
"""
struct ExplicitGridMapping{N,T,G} <: GridMapping{N,T}
    grid::G
    pos2id::Dict{NTuple{N,Int},Int}
    id2pos::Vector{NTuple{N,Int}}
    min_pos::NTuple{N,Int}
    max_pos::NTuple{N,Int}
    n_per_dim::NTuple{N,Int}
    periodic_dims::NTuple{N,Bool}
end

function ExplicitGridMapping(grid, positions::AbstractVector{NTuple{N,Int}};
    periodic_dims::NTuple{N,Bool}=ntuple(_->false,N),
) where {N}

    id2pos = collect(positions)
    pos2id = Dict{NTuple{N,Int},Int}((p,i) for (i,p) in enumerate(id2pos))

    # infer bounds
    mins = ntuple(d -> minimum(p[d] for p in id2pos), N)
    maxs = ntuple(d -> maximum(p[d] for p in id2pos), N)
    n_per_dim = ntuple(d -> maxs[d] - mins[d] + 1, N)

    return ExplicitGridMapping{N, eltype(Dionysos.Domain.get_origin(grid)), typeof(grid)}(
        grid, pos2id, id2pos, mins, maxs, n_per_dim, periodic_dims
    )
end

get_grid(m::ExplicitGridMapping) = m.grid
get_n_state(m::ExplicitGridMapping) = length(m.id2pos)
enum_states(m::ExplicitGridMapping) = 1:get_n_state(m)

is_valid_pos(m::ExplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} =
    haskey(m.pos2id, wrap_pos(m, pos))

get_state_by_pos(m::ExplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} =
    m.pos2id[wrap_pos(m, pos)]
get_pos_by_state(m::ExplicitGridMapping{N}, q::Int) where {N} = m.id2pos[q]

periodic_dims(m::ExplicitGridMapping{N}) where {N} = m.periodic_dims

function wrap_pos(m::ExplicitGridMapping{N}, pos::NTuple{N,Int}) where {N}
    return ntuple(d -> begin
        if !m.periodic_dims[d]
            pos[d]
        else
            span = m.n_per_dim[d]
            z = pos[d] - m.min_pos[d]
            m.min_pos[d] + mod(z, span)
        end
    end, N)
end

"""
ImplicitGridMapping:
- universe is a rectangular index box (min_pos..max_pos) possibly periodic
- id computed by row-major linearization (no dict)
- pos recovered by inverse mapping
"""
struct ImplicitGridMapping{N,T,G,P} <: GridMapping{N,T}
    grid::G
    min_pos::NTuple{N,Int}          # inclusive
    max_pos::NTuple{N,Int}          # inclusive
    n_per_dim::NTuple{N,Int}        # = max_pos - min_pos + 1
    strides::NTuple{N,Int}          # row-major strides
    periodic_dims::NTuple{N,Bool}   # which dims wrap
end

get_grid(m::ImplicitGridMapping) = m.grid

function ImplicitGridMapping(
    grid,
    min_pos::NTuple{N,Int},
    max_pos::NTuple{N,Int};
    periodic_dims::NTuple{N,Bool} = ntuple(_->false, N),
) where {N}
    n_per_dim = ntuple(d -> max_pos[d] - min_pos[d] + 1, N)
    @assert all(n_per_dim .>= 1) "Invalid bounds: max_pos must be >= min_pos"

    # row-major (dim1 fastest)
    strides = ntuple(d -> begin
        s = 1
        for k in 1:(d-1)
            s *= n_per_dim[k]
        end
        s
    end, N)

    return ImplicitGridMapping{N, eltype(Dionysos.Domain.get_origin(grid)), typeof(grid), Nothing}(
        grid, min_pos, max_pos, n_per_dim, strides, periodic_dims
    )
end

get_n_state(m::ImplicitGridMapping) = prod(m.n_per_dim)

enum_states(m::ImplicitGridMapping) = 1:get_n_state(m)

periodic_dims(m::ImplicitGridMapping{N}) where {N} = m.periodic_dims

function wrap_pos(m::ImplicitGridMapping{N}, pos::NTuple{N,Int}) where {N}
    return ntuple(d -> begin
        if !m.periodic_dims[d]
            pos[d]
        else
            # wrap into [min_pos, max_pos]
            span = m.n_per_dim[d]
            # shift to 0..span-1, mod, shift back
            z = pos[d] - m.min_pos[d]
            m.min_pos[d] + mod(z, span)
        end
    end, N)
end

is_valid_pos(m::ImplicitGridMapping{N}, pos::NTuple{N,Int}) where {N} = begin
    p = wrap_pos(m, pos)
    # For periodic dims, always true after wrapping; for non-periodic, must be inside box.
    all(d -> (m.periodic_dims[d] || (m.min_pos[d] <= p[d] <= m.max_pos[d])), 1:N)
end
function get_state_by_pos(m::ImplicitGridMapping{N}, pos::NTuple{N,Int}) where {N}
    p = wrap_pos(m, pos)
    # check non-periodic bounds
    for d in 1:N
        if !m.periodic_dims[d] && !(m.min_pos[d] <= p[d] <= m.max_pos[d])
            throw(DomainError(p, "pos out of implicit universe"))
        end
    end
    # 1 + sum((p[d]-min[d]) * stride[d])
    q = 1
    @inbounds for d in 1:N
        q += (p[d] - m.min_pos[d]) * m.strides[d]
    end
    return q
end

function get_pos_by_state(m::ImplicitGridMapping{N}, q::Int) where {N}
    if !(1 <= q <= get_n_state(m))
        throw(DomainError(q, "state out of implicit universe"))
    end
    z = q - 1
    return ntuple(d -> begin
        span = m.n_per_dim[d]
        # extract digit in mixed radix
        digit = (z ÷ m.strides[d]) % span
        m.min_pos[d] + digit
    end, N)
end

end # module Mappings