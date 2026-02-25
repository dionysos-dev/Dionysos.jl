abstract type AbstractStateSet{N} end

contains_state(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} = error("not implemented")
enum_states(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} = error("not implemented")
get_n_state(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} = length(enum_states(S, m))

add_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} = error("not implemented")
remove_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} = error("not implemented")
empty_states!(S::AbstractStateSet{N}) where {N} = error("not implemented")

function add_states!(S::AbstractStateSet{N}, m::AbstractMapping{N}, states) where {N}
    for q in states
        add_state!(S, m, q)
    end
end

function add_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    add_states!(S, m, states)
    return collect(states)
end

function stateset_from_states(::Type{S}, m::AbstractMapping{N}, states) where {N,S<:AbstractStateSet{N}}
    out = S()
    add_states!(out, m, states)
    return out
end

# Convenient default: ExplicitIdSet (BitSet)
stateset_from_states(m::AbstractMapping{N}, states) where {N} = stateset_from_states(ExplicitIdSet{N}, m, states)

function remove_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    for q in states
        remove_state!(S, m, q)
    end
    return collect(states)
end

@recipe function f(
    tup::Tuple{AbstractStateSet{N}, GridMapping{N,T}};
    dims = [1,2],
    efficient = true,
    label = "",
) where {N,T}
    S, m = tup
    grid = get_grid(m)
    d1 = dims[1]
    d2 = dims[2]

    states = enum_states(S, m)  # iterator
    pos2d = positions2d_from_states(m, states; dims=dims)

    first_series = true
    if !efficient
        # plot each cell (slow)
        for (x,y) in pos2d
            @series begin
                label := first_series ? label : ""
                first_series = false
                return grid, (x,y)
            end
        end
    else
        rects = merge_rectangles_2d(pos2d)
        for r in rects
            @series begin
                label := first_series ? label : ""
                first_series = false
                dims := [1,2]
                return intrect2_to_real_rect(grid, r, d1, d2)
            end
        end
    end
end


# ---------- Helpers ----------
struct UniqueStates{I}
    iter::I
end
function Base.iterate(U::UniqueStates, st=(iterate(U.iter), BitSet()))
    (it, seen) = st
    it === nothing && return nothing
    (val, s2) = it
    while in(val, seen)
        it = iterate(U.iter, s2)
        it === nothing && return nothing
        (val, s2) = it
    end
    push!(seen, val)
    return (val, (iterate(U.iter, s2), seen))
end
unique_states(iter) = UniqueStates(iter)
# ---------- Helpers ----------


mutable struct ExplicitIdSet{N} <: AbstractStateSet{N}
    bits::BitSet
end
ExplicitIdSet{N}() where {N} = ExplicitIdSet{N}(BitSet())

contains_state(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = in(q, S.bits)
enum_states(S::ExplicitIdSet{N}, m::AbstractMapping{N}) where {N} = S.bits
add_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = push!(S.bits, q)
remove_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = delete!(S.bits, q)
empty_states!(S::ExplicitIdSet{N}) where {N} = empty!(S.bits)



# --------------------------


struct MappingSet{N} <: AbstractStateSet{N} end

contains_state(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} = is_valid_state(m, q)
enum_states(::MappingSet{N}, m::AbstractMapping{N}) where {N} = 1:get_n_state(m)
get_n_state(::MappingSet{N}, m::AbstractMapping{N}) where {N} = get_n_state(m)

add_state!(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("MappingSet is read-only")
remove_state!(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("MappingSet is read-only")
empty_states!(::MappingSet{N}) where {N} =
    error("MappingSet is read-only")

# -------------------------- 

mutable struct ImplicitStateSet{N} <: AbstractStateSet{N}
    set::UT.LazySetMinus
end

ImplicitStateSet{N}() where {N} =
    ImplicitStateSet{N}(UT.LazySetMinus(UT.LazyUnion([]), UT.LazyUnion([])))


# --------------------------
# Grid helpers
# --------------------------
_cell_center(m::GridMapping, q::Int) = get_coord_by_state(m, q)
struct CellCornerIter{N,T}
    c::SVector{N,T}
    r::SVector{N,T}
end
Base.iterate(it::CellCornerIter{N,T}, mask::Int=0) where {N,T} = _iterate_corners(it, mask)
function _iterate_corners(it::CellCornerIter{N,T}, mask::Int) where {N,T}
    mask >= (1 << N) && return nothing
    c, r = it.c, it.r
    x = SVector{N,T}(ntuple(i -> ((mask >> (i-1)) & 0x01) == 0 ? c[i]-r[i] : c[i]+r[i], N))
    return (x, mask + 1)
end
function _cell_corner_iter(m::GridMapping, q::Int)
    grid = get_grid(m)
    pos  = get_pos_by_state(m, q)
    c    = get_coord_by_pos(grid, pos)
    r    = get_h(grid) / 2
    return CellCornerIter(c, r)
end

contains_state(S::ImplicitStateSet{N}, m::GridMapping{N}, q::Int) where {N} = contains_state(S, m, q; incl_mode=INNER)
function contains_state(S::ImplicitStateSet{N}, m::GridMapping, q::Int; incl_mode=INNER) where {N}
    if !is_valid_state(m, q)
        return false
    end
    set = S.set
    if incl_mode == CENTER
        c = _cell_center(m, q)
        return UT._point_in_set(set, c)

    elseif incl_mode == INNER
        # conservative: all corners must be in A and not in B
        for xcorner in _cell_corner_iter(m, q)
            UT._point_in_set(set.A, xcorner) || return false
            UT._point_in_set(set.B, xcorner) && return false
        end
        return true

    elseif incl_mode == OUTER
        # sufficient: some sample in A and not in B
        c = _cell_center(m, q)
        if UT._point_in_set(set.A, c) && !UT._point_in_set(set.B, c)
            return true
        end
        for xcorner in _cell_corner_iter(m, q)
            if UT._point_in_set(set.A, xcorner) && !UT._point_in_set(set.B, xcorner)
                return true
            end
        end
        return false

    else
        error("Unknown incl_mode=$incl_mode (expected CENTER/INNER/OUTER)")
    end
end

enum_states(S::ImplicitStateSet{N}, m::AbstractMapping{N}) where {N} = enum_states(S, m ; incl_mode=INNER)
enum_states(S::ImplicitStateSet{N}, m::AbstractMapping{N}; incl_mode=INNER) where {N} = get_states_from_set(m, S.set, incl_mode)

function add_set!(S::ImplicitStateSet{N}, m::AbstractMapping, set) where {N}
    S.set = UT.add_set(S.set, set)
    return
end

function remove_set!(S::ImplicitStateSet{N}, m::AbstractMapping, set) where {N}
    S.set = UT.remove_set(S.set, set)
    return
end

function empty_states!(S::ImplicitStateSet{N}) where {N}
    S.set = UT.LazySetMinus(UT.LazyUnion([]), UT.LazyUnion([]))
end

add_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use add_set!(S, m, rect/union/...)")

remove_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use remove_set!(S, m, rect/union/...)")

# --------------------------

struct UnionStateSet{N,S1<:AbstractStateSet{N},S2<:AbstractStateSet{N}} <: AbstractStateSet{N}
    A::S1
    B::S2
end
contains_state(S::UnionStateSet{N}, m, q::Int) where {N} =
    contains_state(S.A, m, q) || contains_state(S.B, m, q)
enum_states(S::UnionStateSet{N}, m) where {N} =
    unique_states(Iterators.flatten((enum_states(S.A,m), enum_states(S.B,m))))


    
# --------------------------

struct SetMinusStateSet{N,S1<:AbstractStateSet{N},S2<:AbstractStateSet{N}} <: AbstractStateSet{N}
    A::S1
    B::S2
end
contains_state(S::SetMinusStateSet{N}, m, q::Int) where {N} =
    contains_state(S.A, m, q) && !contains_state(S.B, m, q)
enum_states(S::SetMinusStateSet{N}, m) where {N} =
    (q for q in enum_states(S.A,m) if !contains_state(S.B,m,q))


# --------------------------