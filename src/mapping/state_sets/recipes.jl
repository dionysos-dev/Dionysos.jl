# ----------------------------
# Plotting: state set over a grid mapping
# ----------------------------

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
                if proj !== nothing
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
                return intrect2_to_real_rect(grid, r, d1, d2)
            end
        end
    end
end

# ----------------------------
# Plotting helpers
# ----------------------------

struct IntRect2
    lb::NTuple{2, Int}
    ub::NTuple{2, Int}
end

function merge_rectangles_2d(pos2d::AbstractVector{NTuple{2, Int}})
    isempty(pos2d) && return IntRect2[]

    rows = Dict{Int, Vector{Int}}()
    for (x, y) in pos2d
        push!(get!(rows, y, Int[]), x)
    end
    for xs in values(rows)
        sort!(xs)
    end

    intervals = Tuple{Int, Int, Int}[]  # (y, x1, x2)
    ys = sort!(collect(keys(rows)))
    for y in ys
        xs = rows[y]
        i = 1
        while i <= length(xs)
            x1 = xs[i]
            x2 = x1
            while i < length(xs) && xs[i + 1] == x2 + 1
                i += 1
                x2 = xs[i]
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
        y1 = y
        y2 = y
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
        push!(rects, IntRect2((x1, y1), (x2, y2)))
    end
    return rects
end

function intrect2_to_real_rect(grid, r::IntRect2, d1::Int, d2::Int)
    orig = get_origin(grid)
    h = get_h(grid)
    T = eltype(orig)

    lb = SVector{2, T}(
        orig[d1] + r.lb[1]*h[d1] - h[d1]/2,
        orig[d2] + r.lb[2]*h[d2] - h[d2]/2,
    )
    ub = SVector{2, T}(
        orig[d1] + r.ub[1]*h[d1] + h[d1]/2,
        orig[d2] + r.ub[2]*h[d2] + h[d2]/2,
    )
    return UT.box(lb, ub)
end

function project_states_on_dims(
    m::GridMapping{N},
    states;
    dims = [1, 2],
    value_function::Union{Nothing, Function} = nothing,
    reducer::Function = min,
) where {N}
    d1, d2 = dims[1], dims[2]

    if value_function === nothing
        # just dedup by (d1,d2), keep any representative full pos
        seen = Set{NTuple{2, Int}}()
        posN = NTuple{N, Int}[]
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
        proj = Dict{NTuple{2, Int}, Tuple{Float64, NTuple{N, Int}}}()
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
