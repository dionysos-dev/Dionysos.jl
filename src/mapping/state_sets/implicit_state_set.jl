# ----------------------------
# Grid cell geometry helpers (used by ImplicitStateSet membership)
# ----------------------------

_cell_center(m::GridMapping, q::Int) = get_coord_by_state(m, q)

struct CellCornerIter{N, T}
    c::SVector{N, T}
    r::SVector{N, T}
end
Base.iterate(it::CellCornerIter{N, T}, mask::Int = 0) where {N, T} =
    _iterate_corners(it, mask)
function _iterate_corners(it::CellCornerIter{N, T}, mask::Int) where {N, T}
    mask >= (1 << N) && return nothing
    c, r = it.c, it.r
    x = SVector{N, T}(ntuple(i -> ((mask >> (i-1)) & 0x01) == 0 ? c[i]-r[i] : c[i]+r[i], N))
    return (x, mask + 1)
end
function _cell_corner_iter(m::GridMapping, q::Int)
    grid = get_grid(m)
    pos = get_pos_by_state(m, q)
    c = get_coord_by_pos(grid, pos)
    r = get_h(grid) / 2
    return CellCornerIter(c, r)
end

# ----------------------------
# ImplicitStateSet
# ----------------------------

"""
    ImplicitStateSet{N} <: AbstractStateSet{N}

A geometry-backed state set: membership is decided lazily from a stored
`LazySet` and an inclusion mode, rather than from an explicit list of labels.

The `set` field is deliberately abstractly typed (`LazySets.LazySet`): the set is
mutated in place through `add_set!`/`remove_set!`, and each operation may return
a *different* concrete `LazySet` subtype (`EmptySet` → `UnionSetArray` →
`Intersection`, …), so the field cannot be a fixed concrete parameter.
"""
mutable struct ImplicitStateSet{N} <: AbstractStateSet{N}
    set::LazySets.LazySet
    incl_mode::INCL_MODE
end

ImplicitStateSet(set, incl_mode::INCL_MODE) =
    ImplicitStateSet{LazySets.dim(set)}(set, incl_mode)

function ImplicitStateSet(m::AbstractMapping, set, incl_mode::INCL_MODE)
    if is_periodic(m)
        set = UT.set_in_period(
            set,
            get_periodic_dims(m),
            get_periods(m),
            get_periodic_starts(m),
        )
    end
    return ImplicitStateSet{LazySets.dim(set)}(set, incl_mode)
end

ImplicitStateSet{N}() where {N} = ImplicitStateSet{N}(LazySets.EmptySet(N), INNER)

Base.copy(S::ImplicitStateSet{N}) where {N} = ImplicitStateSet{N}(S.set, S.incl_mode)

function contains_state(S::ImplicitStateSet{N}, m::GridMapping{N}, q::Int) where {N}
    if !is_valid_state(m, q)
        return false
    end
    set = S.set
    A = UT.minus_included(set)
    B = UT.minus_hole(set)
    if S.incl_mode == CENTER
        c = _cell_center(m, q)
        return c ∈ set

    elseif S.incl_mode == INNER
        # conservative: every corner must lie in A and outside B
        for xcorner in _cell_corner_iter(m, q)
            xcorner ∈ A || return false
            xcorner ∈ B && return false
        end
        return true

    elseif S.incl_mode == OUTER
        # sufficient: some sample lies in A and outside B
        c = _cell_center(m, q)
        if c ∈ A && c ∉ B
            return true
        end
        for xcorner in _cell_corner_iter(m, q)
            if xcorner ∈ A && xcorner ∉ B
                return true
            end
        end
        return false

    else
        error("Unknown incl_mode=$(S.incl_mode) (expected CENTER/INNER/OUTER)")
    end
end

enum_states(S::ImplicitStateSet{N}, m::AbstractMapping{N}) where {N} =
    get_states_from_set(m, S.set, S.incl_mode)

function add_set!(
    S::ImplicitStateSet{N},
    m::AbstractMapping,
    set,
    incl_mode::INCL_MODE,
) where {N}
    if is_periodic(m)
        set = UT.set_in_period(
            set,
            get_periodic_dims(m),
            get_periods(m),
            get_periodic_starts(m),
        )
    end
    S.set = UT.add_set(S.set, set)
    S.incl_mode = incl_mode
    return
end

function remove_set!(S::ImplicitStateSet{N}, m::AbstractMapping, set) where {N}
    if is_periodic(m)
        set = UT.set_in_period(
            set,
            get_periodic_dims(m),
            get_periods(m),
            get_periodic_starts(m),
        )
    end
    S.set = UT.remove_set(S.set, set)
    return
end

function empty_states!(S::ImplicitStateSet{N}) where {N}
    return S.set = LazySets.EmptySet(N)
end

add_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use add_set!(S, m, rect/union/...)")

remove_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use remove_set!(S, m, rect/union/...)")
