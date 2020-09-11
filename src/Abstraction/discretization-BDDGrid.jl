mutable struct BDDGridDisc{N,T} <: GridDisc{N,T}
    grid::TheGrid{N,T}
    posmin::SVector{N,Int}
end

BDDGridDisc(orig::SVector{N}, h::SVector{N}, posmin::SVector{N,Int}) where N =
    BDDGridDisc(TheGrid(orig, h), posmin)

struct BDDGridCell{N}
    vals::SVector{N,UInt}
end

_celltype(::Type{<:BDDGridDisc{N}}) where N = BDDGridCell{N}

# Returns Int(a + b)
_sum1(a::Int, b::UInt) = b >= -a ? Int(a + b) : -Int(-a - b)
_vals2pos(posmin::SVector{N,Int}, vals::SVector{N,UInt}) where N = _sum1.(posmin, vals)
# Return UInt(a - b) where a >= b
_diff1(a::Int, b::Int) = a < 0 ? UInt(-b) + a : UInt(a) - b
_pos2vals(posmin::SVector{N,Int}, pos::SVector{N,Int}) where N = _diff1.(pos, posmin)
_pos2vals(posmin::SVector{N,Int}, post::NTuple{N,Int}) where N = _diff1.(SVector(post), posmin)
_vals2coord(grid::TheGrid, posmin::SVector{N,Int}, vals::SVector{N,UInt}) where N =
    (pos = _vals2pos(posmin, vals); _grid_pos2coord(grid, pos))
_coord2vals(grid::TheGrid, posmin::SVector{N,Int}, x::SVector) where N =
    (pos = _grid_coord2pos(grid, x); _pos2vals(posmin, pos))

pos2cell(disc::BDDGridDisc) = let posmin = disc.posmin
    pos -> BDDGridCell(_pos2vals(posmin, pos))
end
cell2pos(disc::BDDGridDisc) = let posmin = disc.posmin
    cell -> _vals2pos(posmin, cell.vals)
end

coord2cell(disc::BDDGridDisc) = let grid = disc.grid, posmin = disc.posmin
    x -> BDDGridCell(_coord2vals(grid, posmin, x))
end
cell2coord(disc::BDDGridDisc) = let grid = disc.grid, posmin = disc.posmin
    cell -> _vals2coord(grid, posmin, cell.vals)
end
