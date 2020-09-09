mutable struct ListGridStateSet{S} <: StateSet
    elems::S
end
mutable struct ListGridLabelSet{S} <: LabelSet
    elems::S
end

AddStateSet!(mng::ListGridManager) = ListGridStateSet(Set{statetype(mng)}())
AddLabelSet!(mng::ListGridManager) = ListGridLabelSet(Set{labeltype(mng)}())

# ---
Base.empty!(mng::ListGridManager, symbset::SymbolSet) = (empty!(symbset.elems); symbset)
get_nsymbols(mng::ListGridManager, symbset) = length(symbset.elems)
# Return an iterable of all symbols
enum_symbols(mng::ListGridManager, symbset::SymbolSet) = symbset.elems
get_some_symbol(mng::ListGridManager, symbset) =
    isempty(symbset.elems) ? nothing : first(symbset.elems)
# ---

function create_symbols!(mng::ListGridManager, domain::Domain)
    append!(get_map(mng, _XUType(domain)).coll, domain.cells)
end

add_symbol!(mng::ListGridManager, symbset::SymbolSet, symb) = push!(symbset.elems, symb)

function _add_symbols!(mng::ListGridManager, symbset, domain, XUType)
    get_symbol_try = cell2symbol_try(get_map(mng, XUType))
    elems = symbset.elems
    for cell in domain.cells
        symb, is_valid = get_symbol_try(cell)
        if is_valid
            push!(elems, symb)
        end
    end
end

add_symbols!(mng::Manager, symbset::StateSet, domain::XDomain) =
    _add_symbols!(mng, symbset, domain, XType)
add_symbols!(mng::Manager, symbset::LabelSet, domain::UDomain) =
    _add_symbols!(mng, symbset, domain, UType)

#---
function _add_cells!(mng::ListGridManager, domain, symbset, XUType)
    get_cell = symbol2cell(get_map(mng, XUType))
    cells = domain.cells
    for symb in symbset.elems
        push!(cells, get_cell(symb))
    end
end

add_cells!(mng::Manager, domain::XDomain, symbset::StateSet) =
    _add_cells!(mng, domain, symbset, XType)
add_cells!(mng::Manager, domain::UDomain, symbset::LabelSet) =
    _add_cells!(mng, domain, symbset, UType)
