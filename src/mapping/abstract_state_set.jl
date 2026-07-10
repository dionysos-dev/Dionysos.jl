abstract type AbstractStateSet{N} end

contains_state(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `contains_state` for $(typeof(S))")
enum_states(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} =
    error("implement `enum_states` for $(typeof(S))")
get_n_state(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} =
    length(enum_states(S, m))

add_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `add_state!` for $(typeof(S))")
remove_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `remove_state!` for $(typeof(S))")
empty_states!(S::AbstractStateSet{N}) where {N} =
    error("implement `empty_states!` for $(typeof(S))")

function add_states!(S::AbstractStateSet{N}, m::AbstractMapping{N}, states) where {N}
    for q in states
        add_state!(S, m, q)
    end
end

function add_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    add_states!(S, m, states)
    return collect(states)
end

function stateset_from_states(
    ::Type{S},
    m::AbstractMapping{N},
    states,
) where {N, S <: AbstractStateSet{N}}
    out = S()
    add_states!(out, m, states)
    return out
end

# Convenient default: ExplicitIdSet (BitSet)
stateset_from_states(m::AbstractMapping{N}, states) where {N} =
    stateset_from_states(ExplicitIdSet{N}, m, states)

function remove_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    for q in states
        remove_state!(S, m, q)
    end
    return collect(states)
end

@recipe function f(
    tup::Tuple{<:AbstractStateSet{N}, <:GridMapping{N, T}};
    dims = [1, 2],
    efficient = true,
    label = "",
    value_function = nothing,
    reducer = min,
) where {N, T}
    S, m = tup
    d1, d2 = dims[1], dims[2]

    grid = get_grid(m)

    states = enum_states(S, m)
    proj, posN = project_states_on_dims(
        m,
        states;
        dims = [d1, d2],
        value_function = (value_function isa Function ? value_function : nothing),
        reducer = reducer,
    )

    # If value_function is provided, build a colormap range
    if proj !== nothing
        vals = Float64[v for (v, _) in values(proj)]
        finite_vals = filter(isfinite, vals)
        vmax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
    end

    first_series = true
    if !efficient || value_function !== nothing
        # plot each cell (slow) - needs full N-dim pos
        if proj !== nothing
            posN = sort(posN; by = p -> begin
                key = (p[d1], p[d2])
                v = proj[key][1]
                isfinite(v) ? v : -Inf
            end, rev = true)
        end
        for p in posN
            key = (p[d1], p[d2])
            if proj !== nothing
                v = proj[key][1]
            end

            @series begin
                label := first_series ? label : ""
                first_series = false
                dims := dims
                if proj === nothing
                    # nothing
                else
                    v = proj[key][1]
                    vplot = isfinite(v) ? v : NaN
                    colorbar := true
                    clims := (0.0, vmax)
                    fill_z := vplot
                    seriestype := :shape
                    fillalpha := 1.0
                end

                return grid, p
            end
        end
    else
        pos2d = NTuple{2, Int}[(p[d1], p[d2]) for p in posN]
        rects = merge_rectangles_2d(pos2d)
        for r in rects
            @series begin
                label := first_series ? label : ""
                first_series = false
                dims := [1, 2]
                return intrect2_to_real_rect(grid, r, d1, d2)
            end
        end
    end
end

# ---------- Helpers ----------
# ---------- Helpers ----------
struct UniqueStates{I}
    iter::I
end

Base.IteratorEltype(::Type{UniqueStates{I}}) where {I} = Base.HasEltype()
Base.eltype(::Type{UniqueStates{I}}) where {I} = eltype(I)

Base.IteratorSize(::Type{UniqueStates{I}}) where {I} = Base.SizeUnknown()

function Base.length(U::UniqueStates)
    # O(n) but makes Set(itr) work reliably
    seen = BitSet()
    n = 0
    for v in U
        if !(v in seen)
            push!(seen, v)
            n += 1
        end
    end
    return n
end

function Base.iterate(U::UniqueStates, st = (iterate(U.iter), BitSet()))
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

Base.copy(S::ExplicitIdSet{N}) where {N} = ExplicitIdSet{N}(copy(S.bits))

contains_state(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = in(q, S.bits)
enum_states(S::ExplicitIdSet{N}, m::AbstractMapping{N}) where {N} = S.bits
add_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = push!(S.bits, q)
remove_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    delete!(S.bits, q)
empty_states!(S::ExplicitIdSet{N}) where {N} = empty!(S.bits)

# --------------------------

struct MappingSet{N} <: AbstractStateSet{N} end

Base.copy(::MappingSet{N}) where {N} = MappingSet{N}()

contains_state(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    is_valid_state(m, q)
enum_states(::MappingSet{N}, m::AbstractMapping{N}) where {N} = 1:get_n_state(m)
get_n_state(::MappingSet{N}, m::AbstractMapping{N}) where {N} = get_n_state(m)

add_state!(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("MappingSet is read-only")
remove_state!(::MappingSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("MappingSet is read-only")
empty_states!(::MappingSet{N}) where {N} = error("MappingSet is read-only")

# -------------------------- 

mutable struct ImplicitStateSet{N} <: AbstractStateSet{N}
    set::UT.Region
    incl_mode::INCL_MODE
end

ImplicitStateSet(set, incl_mode::INCL_MODE) =
    ImplicitStateSet{UT.get_dim(set)}(set, incl_mode)

function ImplicitStateSet(m::AbstractMapping, set, incl_mode::INCL_MODE)
    if is_periodic(m)
        set = UT.set_in_period(
            set,
            get_periodic_dims(m),
            get_periods(m),
            get_periodic_starts(m),
        )
    end
    return ImplicitStateSet{UT.get_dim(set)}(set, incl_mode)
end

ImplicitStateSet{N}() where {N} = ImplicitStateSet{N}(UT.empty_region(N), INNER)

Base.copy(S::ImplicitStateSet{N}) where {N} = ImplicitStateSet{N}(S.set, S.incl_mode)

# --------------------------
# Grid helpers
# --------------------------
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

function contains_state(S::ImplicitStateSet{N}, m::GridMapping{N}, q::Int) where {N}
    if !is_valid_state(m, q)
        return false
    end
    set = S.set
    A = UT.minus_included(set)
    B = UT.minus_hole(set)
    if S.incl_mode == CENTER
        c = _cell_center(m, q)
        return UT.point_in_set(set, c)

    elseif S.incl_mode == INNER
        # conservative: every corner must lie in A and outside B
        for xcorner in _cell_corner_iter(m, q)
            UT.point_in_set(A, xcorner) || return false
            UT.point_in_set(B, xcorner) && return false
        end
        return true

    elseif S.incl_mode == OUTER
        # sufficient: some sample lies in A and outside B
        c = _cell_center(m, q)
        if UT.point_in_set(A, c) && !UT.point_in_set(B, c)
            return true
        end
        for xcorner in _cell_corner_iter(m, q)
            if UT.point_in_set(A, xcorner) && !UT.point_in_set(B, xcorner)
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
    return S.set = UT.empty_region(N)
end

add_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use add_set!(S, m, rect/union/...)")

remove_state!(::ImplicitStateSet{N}, m::AbstractMapping, q::Int) where {N} =
    error("ImplicitStateSet is geometry-backed. Use remove_set!(S, m, rect/union/...)")

# --------------------------

struct UnionStateSet{N, S1 <: AbstractStateSet{N}, S2 <: AbstractStateSet{N}} <:
       AbstractStateSet{N}
    A::S1
    B::S2
end

contains_state(S::UnionStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    contains_state(S.A, m, q) || contains_state(S.B, m, q)

enum_states(S::UnionStateSet{N}, m::AbstractMapping{N}) where {N} =
    unique_states(Iterators.flatten((enum_states(S.A, m), enum_states(S.B, m))))

# --------------------------

struct SetMinusStateSet{N, S1 <: AbstractStateSet{N}, S2 <: AbstractStateSet{N}} <:
       AbstractStateSet{N}
    A::S1
    B::S2
end

contains_state(S::SetMinusStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    contains_state(S.A, m, q) && !contains_state(S.B, m, q)

enum_states(S::SetMinusStateSet{N}, m::AbstractMapping{N}) where {N} =
    (q for q in enum_states(S.A, m) if !contains_state(S.B, m, q))
# --------------------------

"""
    MappedStateSet{N, S, M}

A state set bundled with the mapping that gives its integer labels meaning.
Consumers hold one self-contained object instead of threading the
`(set, mapping)` pair through every call — which also removes the
"wrong mapping" bug class.
"""
struct MappedStateSet{N, S <: AbstractStateSet{N}, M <: AbstractMapping{N}}
    set::S
    mapping::M
end

get_set(ms::MappedStateSet) = ms.set
get_mapping(ms::MappedStateSet) = ms.mapping
get_dim(::MappedStateSet{N}) where {N} = N

contains_state(ms::MappedStateSet, q::Int) = contains_state(ms.set, ms.mapping, q)
enum_states(ms::MappedStateSet) = enum_states(ms.set, ms.mapping)
get_n_state(ms::MappedStateSet) = get_n_state(ms.set, ms.mapping)

add_state!(ms::MappedStateSet, q::Int) = add_state!(ms.set, ms.mapping, q)
remove_state!(ms::MappedStateSet, q::Int) = remove_state!(ms.set, ms.mapping, q)
add_states!(ms::MappedStateSet, states) = add_states!(ms.set, ms.mapping, states)
add_set!(ms::MappedStateSet, set, incl_mode::INCL_MODE) =
    add_set!(ms.set, ms.mapping, set, incl_mode)
remove_set!(ms::MappedStateSet, set, incl_mode::INCL_MODE) =
    remove_set!(ms.set, ms.mapping, set, incl_mode)

Base.copy(ms::MappedStateSet) = MappedStateSet(copy(ms.set), ms.mapping)

@recipe function f(ms::MappedStateSet)
    return ((ms.set, ms.mapping),) # delegates to the tuple recipe
end
