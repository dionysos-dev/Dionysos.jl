"""
    AbstractMultiLevelMapping{N, T} <: AbstractMapping{N, T}

A mapping whose labels are distributed across several resolution levels. Unlike a
[`GridMapping`](@ref), a position is a `(level, pos)` pair and the pos-based
methods take an explicit `level` argument — so a multi-level mapping is
deliberately *not* a `GridMapping{N,T}` (it does not satisfy the single-grid
`get_pos_by_state(m, q)::NTuple{N,Int}` contract).
"""
abstract type AbstractMultiLevelMapping{N, T} <: AbstractMapping{N, T} end

"""
    HierarchicalGridMapping{N, T, ML} <: AbstractMultiLevelMapping{N, T}

Multiscale grid mapping with several levels of resolution; each level is itself a
`GridMapping` (e.g. an `ImplicitGridMapping`). Global labels are assigned by
concatenating the levels (`offsets[l]` is the number of labels before level `l`).
"""
mutable struct HierarchicalGridMapping{N, T, ML <: GridMapping{N, T}} <:
               AbstractMultiLevelMapping{N, T}
    levels::Vector{ML}
    offsets::Vector{Int}  # offsets[l] = sum_{k<l} n_k; offsets[1]=0
end

function HierarchicalGridMapping(levels::Vector{ML}) where {N, T, ML <: GridMapping{N, T}}
    h = HierarchicalGridMapping{N, T, ML}(levels, Int[])
    _recompute_offsets!(h)
    return h
end

function HierarchicalGridMapping(
    m0::ImplicitGridMapping{N, T};
    nlevels::Int = 1,
    div::Int = 2,
) where {N, T}
    levels = ImplicitGridMapping{N, T, typeof(get_grid(m0))}[]
    push!(levels, m0)
    for _ in 2:nlevels
        push!(levels, build_refined_level(levels[end]; div = div))
    end
    return HierarchicalGridMapping(levels)
end

function _recompute_offsets!(h::HierarchicalGridMapping)
    L = length(h.levels)
    resize!(h.offsets, L)
    h.offsets[1] = 0
    for l in 2:L
        h.offsets[l] = h.offsets[l - 1] + get_n_state(h.levels[l - 1])
    end
    return h
end

get_n_state(h::HierarchicalGridMapping) = h.offsets[end] + get_n_state(h.levels[end])

# expose the underlying grid at a level
get_grid(h::HierarchicalGridMapping, l::Int) = get_grid(h.levels[l])
get_grid(h::HierarchicalGridMapping) = get_grid(h, 1)

# Locate global id -> (level, local_id)
function _locate(h::HierarchicalGridMapping, q::Int)
    (1 <= q <= get_n_state(h)) || return nothing
    for l in 1:length(h.levels)
        n = get_n_state(h.levels[l])
        if h.offsets[l] < q <= h.offsets[l] + n
            return l, q - h.offsets[l]
        end
    end
    return error("unreachable")  # offsets inconsistent
end

# Global -> coord
function get_coord_by_state(h::HierarchicalGridMapping{N, T}, q::Int) where {N, T}
    l, qloc = _locate(h, q)
    return get_coord_by_state(h.levels[l], qloc)
end

# Global -> elem
function get_elem_by_state(h::HierarchicalGridMapping{N, T}, q::Int) where {N, T}
    l, qloc = _locate(h, q)
    return get_elem_by_state(h.levels[l], qloc)
end

function get_state_by_coord(h::HierarchicalGridMapping{N}, x, l::Int) where {N}
    qloc = get_state_by_coord(h.levels[l], x)
    return h.offsets[l] + qloc
end

function get_state_by_coord(h::HierarchicalGridMapping{N}, x) where {N}
    l = length(h.levels)  # finest
    qloc = get_state_by_coord(h.levels[l], x)
    return h.offsets[l] + qloc
end

# Multi-level position is a (level, pos) pair (not the GridMapping NTuple{N,Int}).
function get_pos_by_state(h::HierarchicalGridMapping{N}, q::Int) where {N}
    l, qloc = _locate(h, q)
    return (l, get_pos_by_state(h.levels[l], qloc))
end

function get_state_by_pos(
    h::HierarchicalGridMapping{N},
    l::Int,
    pos::NTuple{N, Int},
) where {N}
    qloc = get_state_by_pos(h.levels[l], pos)
    return h.offsets[l] + qloc
end

function add_level!(h::HierarchicalGridMapping{N, T, ML}; div::Int = 2) where {N, T, ML}
    last = h.levels[end]
    last isa ImplicitGridMapping || error(
        "add_level! default only supports ImplicitGridMapping; provide a constructor for your mapping type",
    )
    push!(h.levels, build_refined_level(last; div = div))
    push!(h.offsets, h.offsets[end] + get_n_state(last))
    return h
end

function get_states_from_set(h::HierarchicalGridMapping, l::Int, set, incl_mode::INCL_MODE)
    qs = get_states_from_set(h.levels[l], set, incl_mode)
    off = h.offsets[l]
    return [off + q for q in qs]
end

function build_refined_level(m::ImplicitGridMapping{N, T, G}; div::Int = 2) where {N, T, G}
    g = get_grid(m)
    g2 = GridFree(get_origin(g), get_h(g) / div)
    # scale index box to keep same physical coverage
    min2 = ntuple(i -> m.min_pos[i] * div, N)
    max2 = ntuple(i -> (m.max_pos[i] + 1) * div - 1, N)
    return ImplicitGridMapping(g2, min2, max2)
end
