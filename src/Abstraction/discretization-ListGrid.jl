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

pos2cell(disc::ListGridDisc) = pos -> ListGridCell(pos)
cell2pos(disc::ListGridDisc) = cell -> cell.pos

coord2cell(disc::ListGridDisc) = let grid = disc.grid
    x -> ListGridCell(_grid_coord2pos(grid, x))
end
cell2coord(disc::ListGridDisc) = let grid = disc.grid
    cell -> _grid_pos2coord(grid, cell.pos)
end
