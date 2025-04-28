"""
    NestedDomain{N, T}

A multiresolution **hierarchical domain** structure that stores a **stack of grid-based domains** (`GridDomainType`), allowing dynamic **refinement** of the domain in space.

Each level represents a **grid discretization** of the state space, where finer grids can be added progressively.

# Fields
- `domains::Vector{GridDomainType{N, T}}`:  
  List of domain grids at each refinement level.  
- `active::Vector{Dict{NTuple{N, Int}, Bool}}`:  
  Dictionary at each level storing which cells are currently **active** (occupied).
- `levels::Int`:  
  Current number of refinement levels.

# Key Features
- **Dynamic Refinement**:  
  You can subdivide a cell into finer cells on demand (`cut_pos!`).
- **Multi-level Access**:  
  Query which level (`depth`) a point belongs to.
- **Efficient Storage**:  
  Only active cells are stored at each level.
- **Grid-Based**:  
  Works with any domain subtype implementing `GridDomainType`.

# Typical Usage
- Start from a coarse grid domain.
- Add levels by splitting cells when higher resolution is needed.
- Track multi-resolution hierarchical structure easily.
"""
mutable struct NestedDomain{N, T}
    domains::Vector{GridDomainType{N, T}}
    active::Vector{Dict{NTuple{N, Int}, Bool}}
    levels::Int
end

# Constructor
function NestedDomain(dom::GridDomainType{N, T}) where {N, T}
    dict = Dict{NTuple{N, Int}, Bool}()
    for pos in enum_pos(dom)
        dict[pos] = true
    end
    return NestedDomain{N, T}([dom], [dict], 1)
end

# Get number of levels
get_levels(Ndomain::NestedDomain) = Ndomain.levels

# Check if a position is active at level l
function is_active(Ndomain::NestedDomain, pos, l)
    if l > Ndomain.levels
        return false
    end
    return get(Ndomain.active[l], pos, false)
end

# Get position index from coordinate at level l
get_pos_by_coord(Ndomain::NestedDomain, l, x) = get_pos_by_coord(Ndomain.domains[l], x)

# Get coordinate from position at level l
get_coord_by_pos(Ndomain::NestedDomain, l, pos) = get_coord_by_pos(Ndomain.domains[l], pos)

# Check if a pos exists in the domain at level l
is_pos(Ndomain::NestedDomain, pos, l) = pos in Ndomain.domains[l]

# Find which level a point belongs to
function get_depth(Ndomain::NestedDomain, x)
    for l in 1:Ndomain.levels
        pos = get_pos_by_coord(Ndomain, l, x)
        if is_pos(Ndomain, pos, l) && is_active(Ndomain, pos, l)
            return l
        end
    end
    return 0
end

# Add a new domain
function add_dom!(Ndomain::NestedDomain, dom::GridDomainType)
    push!(Ndomain.domains, dom)
    push!(Ndomain.active, Dict{NTuple{length(get_h(dom)), Int}, Bool}())
    return Ndomain.levels += 1
end

# Get the grid at a given level
get_grid(Ndomain::NestedDomain, l::Int) = get_grid(Ndomain.domains[l])

# Get grid for a point x
get_grid(Ndomain::NestedDomain, x::SVector) = get_grid(Ndomain, get_depth(Ndomain, x))

# Add a finer subdomain by halving the current grid step
function add_sub_dom!(Ndomain::NestedDomain{N, T}; div = 2) where {N, T}
    l = Ndomain.levels
    dom = Ndomain.domains[l]
    subdomain = rescale_domain(dom, 1/div)
    push!(Ndomain.domains, subdomain)
    push!(Ndomain.active, Dict{NTuple{N, Int}, Bool}())
    return Ndomain.levels += 1
end

# Given a pos, return the subdivided child positions after cut
function get_subpos(pos, div::Int)
    lbI = Tuple(p * div for p in pos)
    ubI = Tuple(p * div + div - 1 for p in pos)
    rectI = UT.HyperRectangle(lbI, ubI)
    return Iterators.product(_ranges(rectI)...)
end

# Refine a cell by splitting it into subcells
function cut_pos!(Ndomain::NestedDomain, pos, l; div = 2)
    if l == Ndomain.levels
        add_sub_dom!(Ndomain; div = div)
    end
    dict = Ndomain.active[l + 1]
    Ndomain.active[l][pos] = false
    subpositions = get_subpos(pos, div)
    for spos in subpositions
        dict[spos] = true
    end
    return subpositions
end

# Empty all domains
function Base.empty!(Ndomain::NestedDomain)
    for dom in Ndomain.domains
        empty!(dom)
    end
end

# Find pos and level from a coordinate
function get_pos_by_coord(Ndomain::NestedDomain, x)
    l = get_depth(Ndomain, x)
    return (get_pos_by_coord(Ndomain, l, x), l)
end

# Check if the whole NestedDomain is empty
function Base.isempty(Ndomain::NestedDomain)
    return all(isempty, Ndomain.domains)
end

# Number of active positions
function get_ncells(Ndomain::NestedDomain)
    return sum(count -> count[2], Iterators.flatten(Ndomain.active))
end

# Enumerate all positions across all levels
function enum_pos(Ndomain::NestedDomain)
    return flatten([enum_pos(dom) for dom in Ndomain.domains])
end

# Plotting recipe
@recipe function f(nd::NestedDomain)
    for l in 1:get_levels(nd)
        grid = get_grid(nd, l)
        for (pos, v) in nd.active[l]
            if v
                @series begin
                    return grid, pos
                end
            end
        end
    end
end
