abstract type SymbolicModel end

struct SymbolicModelList{NX, NU, S<:Automaton} <: SymbolicModel
    Xgrid::GridSpaceList{NX}
    Ugrid::GridSpaceList{NU}
    automaton::S
    xpos2int::Dict{NTuple{NX, Int}, Int}
    xint2pos::Vector{NTuple{NX, Int}}
    upos2int::Dict{NTuple{NU, Int}, Int}
    uint2pos::Vector{NTuple{NU, Int}}
end

function NewSymbolicModelListList(
    Xgrid::GridSpaceList{NX}, Ugrid::GridSpaceList{NU}) where {NX, NU}
    #---------------------------------------------------------------------------
    nx = get_ncells(Xgrid)
    nu = get_ncells(Ugrid)
    xint2pos = [pos for pos in enum_pos(Xgrid)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Xgrid)))
    uint2pos = [pos for pos in enum_pos(Ugrid)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Ugrid)))
    automaton = NewAutomatonList(nx, nu)
    return SymbolicModelList(Xgrid, Ugrid, automaton, xpos2int, xint2pos, upos2int, uint2pos)
end

function get_xpos_by_state(symmodel::SymbolicModelList, state)
    return symmodel.xint2pos[state]
end

function get_state_by_xpos(symmodel::SymbolicModelList, xpos)
    return symmodel.xpos2int[xpos]
end

function get_upos_by_symbol(symmodel::SymbolicModelList, symbol)
    return symmodel.uint2pos[symbol]
end

function get_symbol_by_upos(symmodel::SymbolicModelList, upos)
    return symmodel.upos2int[upos]
end
