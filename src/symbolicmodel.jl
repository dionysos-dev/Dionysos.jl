function NewSymbolicModelHashList(X_grid::GridSpaceHash, U_grid::GridSpaceHash)
    nx = get_gridspace_size(X_grid)
    nu = get_gridspace_size(U_grid)
    xreflist = [ref for ref in enumerate_gridspace_ref(X_grid)]
    ureflist = [ref for ref in enumerate_gridspace_ref(U_grid)]
    automaton = NewAutomatonList(nx, nu)
    return SymbolicModelHash(X_grid, U_grid, automaton, xreflist, ureflist)
end

function get_xref_by_state(sym_model::SymbolicModelHash, state)
    return sym_model.xreflist[state]
end

function get_uref_by_symbol(sym_model::SymbolicModelHash, symbol)
    return sym_model.ureflist[symbol]
end
