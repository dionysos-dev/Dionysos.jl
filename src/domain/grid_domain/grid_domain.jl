"""
    abstract type GridDomainType{N, T} <: DomainType{N, T}

An abstract interface for grid-based domains, enforcing required methods for interaction.
"""
abstract type GridDomainType{N, T} <: DomainType{N, T} end  # More specific domain type

# ----------------------------
# Required API for Grid Domains
# ----------------------------
function get_grid(domain::GridDomainType) end
function enum_pos(domain::GridDomainType) end
function get_ncells(domain::GridDomainType) end
function get_somepos(domain::GridDomainType, pos) end
function Base.isempty(domain::GridDomainType) end
function Base.in(pos, domain::GridDomainType) end
function Base.issubset(domain1::GridDomainType, domain2::GridDomainType) end

function add_pos!(domain::GridDomainType, pos) end
function Base.union!(domain1::GridDomainType, domain2::GridDomainType) end
function Base.setdiff!(domain1::GridDomainType, domain2::GridDomainType) end
function Base.empty!(domain::GridDomainType) end
function remove_pos!(domain::GridDomainType, pos) end
# Create a new domain based on `domain`, but replacing the grid step size `h` with `new_h`.
function rescale_domain(domain::GridDomainType, scale::Float64) end

# ----------------------------
# Derived Utility Methods
# ----------------------------

get_pos_by_coord(domain::GridDomainType, coord) = get_pos_by_coord(get_grid(domain), coord)
get_coord_by_pos(domain::GridDomainType, pos) = get_coord_by_pos(get_grid(domain), pos)
get_elem_by_pos(domain::GridDomainType, pos) = get_elem_by_pos(get_grid(domain), pos)
get_elem_by_coord(domain::GridDomainType, coord) =
    get_elem_by_pos(domain, get_pos_by_coord(domain, coord))

get_dim(domain::GridDomainType) = get_dim(get_grid(domain))
enum_coords(domain::GridDomainType) =
    [get_coord_by_pos(domain, pos) for pos in enum_pos(domain)]
enum_elems(domain::GridDomainType) =
    [get_elem_by_pos(domain, pos) for pos in enum_pos(domain)]
add_coord!(domain::GridDomainType, x) = add_pos!(domain, get_pos_by_coord(domain, x))
crop_to_domain(domain::GridDomainType, list_pos) = list_pos ∩ enum_pos(domain)
convert_to_custom_domain(domain::GridDomainType) = CustomList(enum_coords(domain))

function get_subset_pos(
    domain::GridDomainType,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(get_grid(domain), rect, incl_mode)
    return [pos for pos in Iterators.product(_ranges(rectI)...) if pos ∈ domain]
end

function get_subset_pos_in_grid(
    domain::GridDomainType,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(get_grid(domain), rect, incl_mode)
    return Iterators.product(_ranges(rectI)...)
end

function add_set!(domain::GridDomainType, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(get_grid(domain), rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        add_pos!(domain, pos)
    end
end

function add_set!(
    domain::GridDomainType,
    unionSetArray::UT.LazyUnionSetArray,
    incl_mode::INCL_MODE,
)
    for set in unionSetArray.sets
        add_set!(domain, set, incl_mode)
    end
end

function add_set!(domain::GridDomainType, setMinus::UT.LazySetMinus, incl_mode::INCL_MODE)
    add_set!(domain, setMinus.A, incl_mode)
    return remove_set!(domain, setMinus.A ∩ setMinus.B, _invInclMode(incl_mode))
end

# add_subset! adds a subset of positions from domain2 to domain1, but only if they fall within the specified hyperrectangle rect
function add_subset!(
    domain1::GridDomainType,
    domain2::GridDomainType,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(get_grid(domain1), rect, incl_mode)
    # Decide whether to iterate over `rectI` or `domain2` for efficiency
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain2)
        for pos in pos_iter
            if pos ∈ domain2
                add_pos!(domain1, pos)
            end
        end
    else
        # Iterate over `domain2` and check if positions are in `rectI`
        for pos in enum_pos(domain2)
            if pos ∈ rectI
                add_pos!(domain1, pos)
            end
        end
    end
end

function remove_coord!(domain::GridDomainType, x)
    return remove_pos!(domain, get_pos_by_coord(domain, x))
end

function remove_set!(
    domain::GridDomainType,
    unionSetArray::UT.LazyUnionSetArray,
    incl_mode::INCL_MODE,
)
    for set in unionSetArray.sets
        remove_set!(domain, set, incl_mode)
    end
end

function remove_set!(domain::GridDomainType, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    posL = get_subset_pos(domain, rect, incl_mode)
    if length(posL) < get_ncells(domain)
        for pos in posL
            remove_pos!(domain, pos)
        end
    else
        for pos in enum_pos(domain)
            if pos ∈ posL
                remove_pos!(domain, pos)
            end
        end
    end
end

function merge_rectangles_2d(pos2d::AbstractVector{NTuple{2, Int}})
    isempty(pos2d) && return UT.HyperRectangle{NTuple{2, Int}}[]

    # Group by y (row) and collect x's
    rows = Dict{Int, Vector{Int}}()
    for (x, y) in pos2d
        push!(get!(rows, y, Int[]), x)
    end
    for xs in values(rows)
        sort!(xs)
    end

    # Build horizontal intervals per row: (y, x1, x2)
    intervals = Vector{Tuple{Int, Int, Int}}()
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

    # Merge identical intervals vertically
    sort!(intervals)  # sorts by y then x1 then x2
    rects = UT.HyperRectangle{NTuple{2, Int}}[]
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
        push!(rects, UT.HyperRectangle((x1, y1), (x2, y2)))
    end

    return rects
end

function merge_to_hyperrectangles_pos_2d(domain::GridDomainType; dims = [1, 2])
    pos2d = NTuple{2, Int}[]
    seen = Set{NTuple{2, Int}}()
    for p in enum_pos(domain)
        q = (p[dims[1]], p[dims[2]])
        if !(q in seen)
            push!(seen, q)
            push!(pos2d, q)
        end
    end
    return merge_rectangles_2d(pos2d)
end

function merge_to_hyperrectangles_real_2d(domain::GridDomainType; dims = [1, 2])
    grid = get_grid(domain)
    rects_pos = merge_to_hyperrectangles_pos_2d(domain; dims = dims)

    orig = get_origin(grid)
    h = get_h(grid)
    d1 = dims[1]
    d2 = dims[2]
    T = eltype(orig)

    return [
        UT.HyperRectangle(
            SVector{2, T}(
                orig[d1] + r.lb[1]*h[d1] - h[d1]/2,
                orig[d2] + r.lb[2]*h[d2] - h[d2]/2,
            ),
            SVector{2, T}(
                orig[d1] + r.ub[1]*h[d1] + h[d1]/2,
                orig[d2] + r.ub[2]*h[d2] + h[d2]/2,
            ),
        ) for r in rects_pos
    ]
end

@recipe function f(
    domain::GridDomainType{N, T};
    dims = [1, 2],
    efficient = true,
    label = "",
) where {N, T}
    first_series = true
    if !efficient
        grid = get_grid(domain)
        dict = Dict{NTuple{2, Int}, Any}()
        for pos in enum_pos(domain)
            if !haskey(dict, pos[dims])
                dict[pos[dims]] = true
                @series begin
                    label := first_series ? label : ""
                    first_series = false
                    dims := dims
                    return grid, pos
                end
            end
        end
    else
        for real_hyperrect in merge_to_hyperrectangles_real_2d(domain; dims = dims)
            @series begin
                label := first_series ? label : ""
                first_series = false
                dims := [1, 2]   # IMPORTANT: rect is already 2D
                return real_hyperrect
            end
        end
    end
end
