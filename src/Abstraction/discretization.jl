abstract type GridDisc{N,T} <: Discretization end

_spacetype(::Type{<:GridDisc{N,T}}) where {N,T} = (N, T)

# Simple struct representing a uniform grid
struct TheGrid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
end

mutable struct ListGridDisc{N,T} <: GridDisc{N,T}
    grid::TheGrid{N,T}
end

ListGridDisc(orig::VT, h::VT) where {VT<:SVector} = ListGridDisc(TheGrid(orig, h))

struct ListGridCell{N}
    pos::SVector{N,Int}
end

const ListGridCell_seed = hash("ListGridCell")
@generated function Base.hash(cell::ListGridCell{N}, h::UInt) where N
    quote
        h += ListGridCell_seed
        @nexprs $N i -> h = hash(cell.pos[i], h)
    end
end

ListGridCell(post::NTuple{N,Int}) where N = ListGridCell(SVector(post))

_celltype(::Type{<:ListGridDisc{N}}) where N = ListGridCell{N}

_round_trunc(round_mode, x::T) where T = x >= T(INT_MAX) ? INT_MAX :
    x < T(INT_MIN) ? INT_MIN : round(Int, x, round_mode)
_ceil_trunc(x) = _round_trunc(RoundUp, x)
_floor_trunc(x) = _round_trunc(RoundDown, x)
_round_trunc(x) = _round_trunc(RoundNearest, x)

# compute (max(a - b, 0), safe) with safe = true if a - b <= INT_MAX
function _diff_trunc(a::Int, b::Int)
    a <= b && return (0, true)
    # Here a > b
    (a < 0 || b > 0) && return (a - b, true)
    d = a - b
    return (d, d >= 0)
end

_grid_pos2coord(grid::TheGrid, pos::SVector{N,Int}) where N = grid.orig + pos.*grid.h
_grid_coord2pos(grid::TheGrid, x::SVector) = round.(Int, (x - grid.orig)./grid.h)
function _grid_coord2pos_set(grid::TheGrid,
        rect::HyperRectangle{<:SVector}, incl_mode::INCL_MODE)
    α::Float64 = incl_mode == INNER ? 0.5 : -0.5
    lbI = _ceil_trunc.((rect.lb - grid.orig)./grid.h .+ α)
    ubI = _floor_trunc.((rect.ub - grid.orig)./grid.h .- α)
    return HyperRectangle(lbI, ubI)
end

coord2pos(disc::ListGridDisc) = let grid = disc.grid
    x -> _grid_coord2pos(grid, x)
end
coord2pos_set(disc::ListGridDisc) = let grid = disc.grid
    (set, incl_mode) -> _grid_coord2pos_set(grid, set, incl_mode)
end
pos2coord(disc::ListGridDisc) = let grid = disc.grid
    pos -> _grid_pos2coord(grid, pos)
end

pos2cell(disc::ListGridDisc) = pos -> ListGridCell(pos)
cell2pos(disc::ListGridDisc) = cell -> cell.pos

coord2cell(disc::ListGridDisc) = let grid = disc.grid
    x -> ListGridCell(_grid_coord2pos(grid, x))
end
cell2coord(disc::ListGridDisc) = let grid = disc.grid
    cell -> _grid_pos2coord(grid, cell.pos)
end
