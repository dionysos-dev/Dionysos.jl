abstract type AbstractStateSet{N} end

contains_state(S::AbstractStateSet{N}, m, q::Int) where {N} = error("not implemented")
enum_states(S::AbstractStateSet{N}, m) where {N} = error("not implemented")

add_state!(S::AbstractStateSet{N}, m, q::Int) where {N} = error("not implemented")
remove_state!(S::AbstractStateSet{N}, m, q::Int) where {N} = error("not implemented")
empty_states!(S::AbstractStateSet{N}) where {N} = error("not implemented")

function add_set!(
    S::AbstractStateSet{N},
    m,
    set,
    incl_mode::DO.INCL_MODE;
    return_states::Bool = false,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    for q in states
        add_state!(S, m, q)
    end
    return return_states ? collect(states) : S
end

function remove_set!(
    S::AbstractStateSet{N},
    m,
    set,
    incl_mode::DO.INCL_MODE;
    return_states::Bool = false,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    for q in states
        remove_state!(S, m, q)
    end
    return return_states ? collect(states) : S
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

contains_state(S::ExplicitIdSet{N}, m, q::Int) where {N} = in(q, S.bits)
enum_states(S::ExplicitIdSet{N}, m) where {N} = S.bits
add_state!(S::ExplicitIdSet{N}, m, q::Int) where {N} = (push!(S.bits, q); S)
remove_state!(S::ExplicitIdSet{N}, m, q::Int) where {N} = (delete!(S.bits, q); S)
empty_states!(S::ExplicitIdSet{N}) where {N} = (empty!(S.bits); S)


# --------------------------


struct MappingSet{N} <: AbstractStateSet{N} end

contains_state(::MappingSet{N}, m, q::Int) where {N} = is_valid_state(m, q)
enum_states(::MappingSet{N}, m) where {N} = 1:get_n_state(m)

add_state!(::MappingSet{N}, m, q::Int) where {N} =
    error("MappingSet is read-only")
remove_state!(::MappingSet{N}, m, q::Int) where {N} =
    error("MappingSet is read-only")
empty_states!(::MappingSet{N}) where {N} =
    error("MappingSet is read-only")

# -------------------------- 

mutable struct ImplicitStateSet{N} <: AbstractStateSet{N}
    sm::UT.LazySetMinus   # sm.A is "added" union, sm.B is "removed" union
end

ImplicitStateSet{N}() where {N} =
    ImplicitStateSet{N}(UT.LazySetMinus(UT.LazyUnionSetArray([]), UT.LazyUnionSetArray([])))

# --------------------------
# Internal helpers: push sets into LazyUnionSetArray
# --------------------------
function _push_union!(U::UT.LazyUnionSetArray, s)
    push!(U.sets, s)
    return U
end

_point_in_set(set, x) = (x ∈ set)
function _point_in_set(unionS::UT.LazyUnionSetArray, x)
    for s in unionS.sets
        _point_in_set(s, x) && return true
    end
    return false
end
function _point_in_set(setminus::UT.LazySetMinus, x)
    return _point_in_set(setminus.A, x) && !_point_in_set(setminus.B, x)
end

# --------------------------
# Grid helpers
# --------------------------
_cell_center(m, q::Int) = DO.get_coord_by_pos(get_grid(m), get_pos_by_state(m, q))
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
function _cell_corner_iter(m, q::Int)
    grid = get_grid(m)
    pos  = get_pos_by_state(m, q)
    c    = DO.get_coord_by_pos(grid, pos)
    r    = DO.get_h(grid) / 2
    return CellCornerIter(c, r)
end

function contains_state(S::ImplicitStateSet{N}, m, q::Int; incl_mode=DO.CENTER) where {N}
    if !is_valid_state(m, q)
        return false
    end
    sm = S.sm
    if incl_mode == DO.CENTER
        c = _cell_center(m, q)
        return _point_in_set(sm, c)

    elseif incl_mode == DO.INNER
        # conservative: all corners must be in A and not in B
        for xcorner in _cell_corner_iter(m, q)
            _point_in_set(sm.A, xcorner) || return false
            _point_in_set(sm.B, xcorner) && return false
        end
        return true

    elseif incl_mode == DO.OUTER
        # sufficient: some sample in A and not in B
        c = _cell_center(m, q)
        if _point_in_set(sm.A, c) && !_point_in_set(sm.B, c)
            return true
        end
        for xcorner in _cell_corner_iter(m, q)
            if _point_in_set(sm.A, xcorner) && !_point_in_set(sm.B, xcorner)
                return true
            end
        end
        return false

    else
        error("Unknown incl_mode=$incl_mode (expected DO.CENTER/DO.INNER/DO.OUTER)")
    end
end

enum_states(::ImplicitStateSet{N}, m) where {N} =
    error("ImplicitStateSet cannot enumerate without a bounding set; use enum_states(S, m, bounding)")

function enum_states(S::ImplicitStateSet{N}, m, bounding::AbstractStateSet{N}; incl_mode=DO.CENTER) where {N}
    return (q for q in enum_states(bounding, m) if contains_state(S, m, q; incl_mode=incl_mode))
end

function add_set!(S::ImplicitStateSet{N}, m, rect::UT.HyperRectangle, incl_mode::DO.INCL_MODE) where {N}
    _push_union!(S.sm.A, rect)
    return S
end

function add_set!(S::ImplicitStateSet{N}, m, unionS::UT.LazyUnionSetArray, incl_mode::DO.INCL_MODE) where {N}
    for s in unionS.sets
        _push_union!(S.sm.A, s)
    end
    return S
end

function add_set!(S::ImplicitStateSet{N}, m, setminus::UT.LazySetMinus, incl_mode::DO.INCL_MODE) where {N}
    _push_union!(S.sm.A, setminus)
    return S
end

function remove_set!(S::ImplicitStateSet{N}, m, rect::UT.HyperRectangle, incl_mode::DO.INCL_MODE) where {N}
    _push_union!(S.sm.B, rect)
    return S
end

function remove_set!(S::ImplicitStateSet{N}, m, unionS::UT.LazyUnionSetArray, incl_mode::DO.INCL_MODE) where {N}
    for s in unionS.sets
        _push_union!(S.sm.B, s)
    end
    return S
end

function empty_states!(S::ImplicitStateSet{N}) where {N}
    S.sm = UT.LazySetMinus(UT.LazyUnionSetArray([]), UT.LazyUnionSetArray([]))
    return S
end

add_state!(::ImplicitStateSet{N}, m, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use add_set!(S, m, rect/union/...)")

remove_state!(::ImplicitStateSet{N}, m, q::Int) where {N} =
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