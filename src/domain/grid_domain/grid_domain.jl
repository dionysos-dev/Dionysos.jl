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

# ----------------------------
# Derived Utility Methods
# ----------------------------

get_pos_by_coord(domain::GridDomainType, coord) = get_pos_by_coord(get_grid(domain), coord)
get_coord_by_pos(domain::GridDomainType, pos) = get_coord_by_pos(get_grid(domain), pos)
get_elem_by_pos(domain::GridDomainType, pos) = get_elem_by_pos(get_grid(domain), pos)

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
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain)
        for pos in pos_iter
            remove_pos!(domain, pos)
        end
    else
        for pos in enum_pos(domain)
            if pos ∈ rectI
                remove_pos!(domain, pos)
            end
        end
    end
end

"""
    merge_to_hyperrectangles_pos(domain::GridDomainType)

Aggregates the grid elements in `GridDomainType` into the largest possible hyperrectangles.
Returns a list of `HyperRectangle`s that efficiently represent the domain.
"""
function merge_to_hyperrectangles_pos(domain::GridDomainType)
    elements = collect(enum_pos(domain))  # Convert set to sorted list
    sort!(elements)

    N = get_dim(domain)  # Retrieve the number of dimensions

    hyperrectangles = Vector{UT.HyperRectangle{NTuple{N, Int}}}()
    visited = Dict(e => false for e in elements)  # Dict to track visited elements

    for e in elements
        if visited[e]
            continue  # Skip already processed elements
        end

        # Initialize hyperrectangle bounds
        lower_bound = Tuple(e)
        upper_bound = Tuple(e)

        # Try to grow the hyperrectangle along each dimension
        for dim in 1:N
            while true
                next_upper = ntuple(i -> upper_bound[i] + (i == dim), N)

                if haskey(visited, next_upper) &&
                   !visited[next_upper] &&
                   all(
                       haskey(visited, ntuple(i -> next_upper[i] - (i == dim), N)) for
                       i in 1:N
                   )
                    upper_bound = next_upper
                    visited[next_upper] = true
                else
                    break
                end
            end
        end

        visited[e] = true
        push!(hyperrectangles, UT.HyperRectangle(lower_bound, upper_bound))
    end
    return hyperrectangles
end

"""
    merge_hyperrectangles_real(domain_list::GridDomainType)

Uses the existing `merge_to_hyperrectangles(domain_list)` function to merge elements
in grid positions, then converts the result into real-world `HyperRectangle`s using `get_rec(grid, pos)`.
"""
function merge_to_hyperrectangles_real(domain::GridDomainType)
    grid = get_grid(domain)

    # Step 1: Get merged hyperrectangles in terms of grid positions
    merged_pos_rects = merge_to_hyperrectangles_pos(domain)

    # Step 2: Convert grid-position-based hyperrectangles to real-world ones
    real_hyperrects = [
        UT.HyperRectangle(get_rec(grid, rect.lb).lb, get_rec(grid, rect.ub).ub) for
        rect in merged_pos_rects
    ]
    return real_hyperrects
end

# grid based domain
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
                    return grid, pos
                end
            end
        end
    else
        for real_hyperrect in merge_to_hyperrectangles_real(domain)
            @series begin
                label := first_series ? label : ""
                first_series = false
                return real_hyperrect
            end
        end
    end
end
