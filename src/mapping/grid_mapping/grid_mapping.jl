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
is_state_cover(m::GridMapping) = is_state_cover(get_grid(m))

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
