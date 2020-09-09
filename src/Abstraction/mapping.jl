mutable struct ListGridMap{CT} <: Mapping
    coll::IndColl{CT}
end

struct ListGridSymbol
    val::UInt
end

const ListGridSymbol_seed = hash("ListGridSymbol")
Base.hash(symb::ListGridSymbol, h::UInt) = hash(symb.val, h + ListGridSymbol_seed)

ListGridMap(::Type{CT}) where CT = ListGridMap(IndColl{CT}())

_symboltype(::Type{<:ListGridMap}) = ListGridSymbol
_nsymbols(map::ListGridMap) = length(map.coll)
# Return an iterable of all symbols
enum_symbols(mng::ListGridManager, map::Mapping) = map.coll.i2e
get_some_symbol(mng::ListGridManager, map::Mapping) =
    isempty(map.coll.i2e) ? nothing : ListGridSymbol(UInt(1))

cell2symbol(map::ListGridMap) = let e2i = map.coll.e2i
    cell -> ListGridSymbol(e2i[cell])
end

function _cell2symbol_try(e2i, cell)
    val = get(e2i, cell, UInt(0))
    return (ListGridSymbol(val), val > 0)
end
cell2symbol_try(map::ListGridMap) = let e2i = map.coll.e2i
    cell -> _cell2symbol_try(e2i, cell)
end

symbol2cell(map::ListGridMap) = let i2e = map.coll.i2e
    symb -> i2e[symb.val]
end
