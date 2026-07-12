"""
GridMapping{N,T} <: AbstractMapping{N,T}:
- adds grid, and pos <-> state conversions
"""
abstract type GridMapping{N, T} <: AbstractMapping{N, T} end

get_grid(m::GridMapping) = error("implement `get_grid` for $(typeof(m))")
get_state_by_pos(m::GridMapping{N}, pos::NTuple{N, Int}) where {N} =
    error("implement `get_state_by_pos` for $(typeof(m))")
get_pos_by_state(m::GridMapping{N}, q::Int) where {N} =
    error("implement `get_pos_by_state` for $(typeof(m))")
is_valid_pos(m::GridMapping{N}, pos::NTuple{N, Int}) where {N} = true

enum_pos(m::GridMapping{N, T}) where {N, T} =
    (get_pos_by_state(m, q) for q in enum_states(m))
get_elem_by_pos(m::GridMapping{N, T}, pos::NTuple{N, Int}) where {N, T} =
    get_elem_by_pos(get_grid(m), pos)

get_elem_by_state(m::GridMapping{N, T}, q::Int) where {N, T} =
    get_elem_by_pos(m, get_pos_by_state(m, q))
get_elem_by_coord(m::GridMapping{N, T}, x) where {N, T} =
    get_elem_by_pos(m, get_pos_by_coord(get_grid(m), x))
get_state_by_coord(m::GridMapping, x) =
    get_state_by_pos(m, get_pos_by_coord(get_grid(m), x))
get_coord_by_state(m::GridMapping{N, T}, q::Int) where {N, T} =
    get_coord_by_pos(get_grid(m), get_pos_by_state(m, q))

function get_states_by_coord(m::GridMapping{N}, x) where {N}
    grid = get_grid(m)
    qs = Int[]
    for pos in get_all_pos_by_coord(grid, x)
        is_valid_pos(m, pos) || continue
        q = get_state_by_pos(m, pos)
        push!(qs, q)
    end
    unique!(qs)
    return qs
end
has_overlapping_cells(m::GridMapping) = has_overlapping_cells(get_grid(m))

convert_to_list_mapping(m::GridMapping) = ListMapping(collect(enum_coords(m)))

# returns (states, allin) where allin=false if any covered grid-pos is invalid
function get_states_from_set_strict(
    m::AbstractMapping{N, T},
    x::AbstractVector,
    ::INCL_MODE,
) where {N, T}
    q = get_state_by_coord(m, x)
    return q===nothing ? (nothing, false) : (Int[q], true)
end

# Any bounded `LazySet` works here: `get_pos_from_set` enumerates the covered
# cells (exact index ranges for boxes, bounding-box candidates certified per
# inclusion mode otherwise).
function get_states_from_set_strict(
    m::GridMapping{N},
    S::LazySets.LazySet,
    incl_mode::INCL_MODE,
) where {N}
    grid = get_grid(m)

    qs = Int[]
    allin = true

    for pos in get_pos_from_set(grid, S, incl_mode)
        p = pos::NTuple{N, Int}
        if is_valid_pos(m, p)
            push!(qs, get_state_by_pos(m, p))
        else
            allin = false
        end
    end
    return qs, allin
end

get_states_from_set_strict(m::GridMapping, ::LazySets.EmptySet, incl_mode::INCL_MODE) =
    (Int[], true)

function get_states_from_set_strict(
    m::GridMapping{N},
    subsets::LazySets.UnionSetArray,
    incl_mode::INCL_MODE,
) where {N}
    acc = Int[]
    allin = true
    for subset in subsets.array
        qs, ok = get_states_from_set_strict(m, subset, incl_mode)
        Base.append!(acc, qs)
        allin &= ok
    end
    unique!(acc)
    return acc, allin
end

function get_states_from_set_strict(
    m::GridMapping{N},
    set::UT.SetMinus,
    incl_mode::INCL_MODE,
) where {N}
    states_A, okA = get_states_from_set_strict(m, UT.minus_included(set), incl_mode)
    states_B, okB =
        get_states_from_set_strict(m, UT.minus_hole(set), invert_incl_mode(incl_mode))
    return setdiff(states_A, states_B), (okA && okB)
end

# The non-strict form just drops the `allin` flag; every set kind (vector,
# LazySet, EmptySet, UnionSetArray, SetMinus) routes through the `_strict`
# variant above, so there is a single generic method here rather than one per
# set type.
get_states_from_set(m::AbstractMapping, set, incl_mode::INCL_MODE) =
    first(get_states_from_set_strict(m, set, incl_mode))
