# Cell-aligned regions and predicate-based domain carving.
#
# Carving a domain (removing the cells where some constraint may be violated)
# used to be expressed geometrically: build slightly-shrunk copies of the
# removed cell boxes and subtract their union, the ε-shrink being a workaround
# for the closed-boundary OUTER discretization swallowing face-adjacent
# neighbors. `CellUnion` removes the workaround at the root: it *is* a set of
# whole cells of a given grid, and its discretization is resolved by dispatch —
# exactly the member cells, under every inclusion mode.

"""
    CellUnion(grid, positions)

The union of whole cells of `grid` at the given `positions` (an iterable of
integer tuples), as a `LazySets.LazySet`.

Its grid discretization is exact by dispatch: under **every** inclusion mode,
`get_pos_from_set` returns exactly the member positions — shared faces are
owned by the member cell, so a cell-aligned hole removes exactly its cells
(no face-adjacent over-removal) and a cell-aligned target is recovered exactly
under `INNER`. Coordinate membership (`x ∈ S`) quantizes `x` to its cell;
position membership (`pos ∈ S`) tests the tuple directly.

Use [`cells_where`](@ref) to build one from a predicate and
[`image_blocked_cells`](@ref) for constraints on a nonlinear image of the
state.
"""
struct CellUnion{N, T, G <: Grid{N, T}} <: LazySets.LazySet{T}
    grid::G
    positions::Set{NTuple{N, Int}}
end

CellUnion(grid::Grid{N}, positions) where {N} =
    CellUnion(grid, Set{NTuple{N, Int}}(positions))

LazySets.dim(::CellUnion{N}) where {N} = N

# LazySets' generic `copy` reflects over fields and would require `copy` on the
# (immutable) grid; the grid is shared, only the position set is owned.
Base.copy(S::CellUnion) = CellUnion(S.grid, copy(S.positions))

Base.isempty(S::CellUnion) = isempty(S.positions)
Base.length(S::CellUnion) = length(S.positions)
Base.iterate(S::CellUnion, args...) = iterate(S.positions, args...)
Base.eltype(::Type{<:CellUnion{N}}) where {N} = NTuple{N, Int}

Base.in(pos::NTuple{N, Int}, S::CellUnion{N}) where {N} = pos in S.positions
Base.in(x::AbstractVector, S::CellUnion) = get_pos_by_coord(S.grid, x) in S.positions

UT._is_empty_region(S::CellUnion) = isempty(S.positions)

function UT._outer_box(S::CellUnion{N}) where {N}
    isempty(S) && error("cannot compute the outer box of an empty CellUnion")
    lo = reduce((a, b) -> min.(a, b), S.positions)
    hi = reduce((a, b) -> max.(a, b), S.positions)
    h = get_h(S.grid)
    return LazySets.Hyperrectangle(;
        low = get_coord_by_pos(S.grid, lo) .- h ./ 2,
        high = get_coord_by_pos(S.grid, hi) .+ h ./ 2,
    )
end

_same_grid(a::Grid, b::Grid) =
    a === b || (get_origin(a) == get_origin(b) && get_h(a) == get_h(b))

function get_pos_from_set(grid::Grid{N}, S::CellUnion{N}, ::INCL_MODE) where {N}
    _same_grid(grid, S.grid) || error(
        "this CellUnion is aligned with a different grid; re-express it on the " *
        "discretization grid before discretizing.",
    )
    return S.positions
end

"""
    cells_where(pred, grid, region, incl_mode = INNER)

The [`CellUnion`](@ref) of the cells of `region` (discretized with
`incl_mode`) whose position satisfies `pred(pos)`. Supports `do`-block syntax;
convert a position to its center with `get_coord_by_pos(grid, pos)`:

```julia
blocked = MP.cells_where(grid, X) do pos
    is_infeasible(MP.get_coord_by_pos(grid, pos))
end
```
"""
function cells_where(pred, grid::Grid{N}, region, incl_mode::INCL_MODE = INNER) where {N}
    return CellUnion(
        grid,
        Set{NTuple{N, Int}}(
            pos for pos in get_pos_from_set(grid, region, incl_mode) if pred(pos)
        ),
    )
end

"""
    image_blocked_cells(g, gradient_bound, grid, region, obstacles, incl_mode = INNER)

Sound pullback of constraints living on a nonlinear image of the state: the
[`CellUnion`](@ref) of every cell of `region` in which *some* point **might**
map into one of the `obstacles` (a `LazySet` or a vector of them) under `g`.

The test is center-based and certified by a Lipschitz bound: with
`|∂gⱼ/∂xᵢ| ≤ gradient_bound[i]`, the image of any point of a cell stays within
`dev = Σᵢ gradient_bound[i] · hᵢ/2` (per component) of the image of its
center, so a cell whose center image clears every obstacle by `dev` is
certified feasible — a *kept* cell never maps into an obstacle. Conservative
(the deviation ball is over-approximated by a box).
"""
function image_blocked_cells(
    g,
    gradient_bound::AbstractVector,
    grid::Grid,
    region,
    obstacles,
    incl_mode::INCL_MODE = INNER,
)
    obstacle_list = obstacles isa AbstractVector ? obstacles : [obstacles]
    boxes = [UT._outer_box(ob) for ob in obstacle_list]
    dev = sum(gradient_bound .* (get_h(grid) ./ 2))

    return cells_where(grid, region, incl_mode) do pos
        y = g(get_coord_by_pos(grid, pos))
        ball = LazySets.Hyperrectangle(y, SVector{length(y)}(fill(dev, length(y))))
        return any(zip(obstacle_list, boxes)) do (ob, box)
            near = all(
                abs(y[i] - LazySets.center(box)[i]) <= dev + box.radius[i] for
                i in eachindex(y)
            )
            return near && !UT.is_disjoint(ball, ob)
        end
    end
end
