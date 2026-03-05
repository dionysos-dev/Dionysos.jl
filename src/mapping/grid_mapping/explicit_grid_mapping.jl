"""
ExplicitGridMapping:
explicit enumeration of positions -> ids (dict) and ids -> positions (vector)
"""
struct ExplicitGridMapping{N, T, G} <: GridMapping{N, T}
    grid::G
    pos2id::Dict{NTuple{N, Int}, Int}
    id2pos::Vector{NTuple{N, Int}}
end

function _explicit_mapping_ctor(::Val{N}, ::Type{T}, grid::G) where {N, T, G}
    return ExplicitGridMapping{N, T, G}(grid, Dict{NTuple{N, Int}, Int}(), NTuple{N, Int}[])
end

# Infer N,T from the grid type
function ExplicitGridMapping(grid::Grid{N, T}) where {N, T}
    return _explicit_mapping_ctor(Val(N), T, grid)
end

function ExplicitGridMapping{N, T}(grid::G, positions) where {N, T, G}
    m = ExplicitGridMapping{N, T}(grid)
    for pos in positions
        add_pos!(m, pos)   # dedup via pos2id
    end
    return m
end

function ExplicitGridMapping{N, T}(grid::G, set, incl_mode::INCL_MODE) where {N, T, G}
    m = ExplicitGridMapping{N, T}(grid)
    add_set!(m, set, incl_mode)
    return m
end

function (::Type{ExplicitGridMapping{N, T}})(grid::G) where {N, T, G}
    return ExplicitGridMapping{N, T, G}(grid, Dict{NTuple{N, Int}, Int}(), NTuple{N, Int}[])
end

get_grid(m::ExplicitGridMapping) = m.grid
get_n_state(m::ExplicitGridMapping) = length(m.id2pos)
enum_states(m::ExplicitGridMapping) = 1:get_n_state(m)

is_valid_pos(m::ExplicitGridMapping{N}, pos::NTuple{N, Int}) where {N} =
    haskey(m.pos2id, pos)
get_state_by_pos(m::ExplicitGridMapping{N}, pos::NTuple{N, Int}) where {N} =
    is_valid_pos(m, pos) ? m.pos2id[pos] : nothing
get_pos_by_state(m::ExplicitGridMapping{N}, q::Int) where {N} = m.id2pos[q]

function add_pos!(m::ExplicitGridMapping{N, T}, pos::NTuple{N, Int}) where {N, T}
    if haskey(m.pos2id, pos)
        return m.pos2id[pos]
    end
    q = length(m.id2pos) + 1
    push!(m.id2pos, pos)
    m.pos2id[pos] = q
    return q
end

function add_set!(m::ExplicitGridMapping{N, T}, set, incl_mode::INCL_MODE) where {N, T}
    grid = get_grid(m)
    for pos in get_pos_from_set(grid, set, incl_mode)
        add_pos!(m, pos)   # dedup happens here
    end
    return m
end
