
function LazySymbolicModel(Xdom::DO.GeneralDomainList{N,DO.RectanglularObstacles{NTuple{N,T1}}}, Udom::DO.DomainType{N,T2}) where {N,M,T1,T2}
    nu = DO.get_ncells(Udom)
    uint2pos = [pos for pos in DO.enum_pos(Udom)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Udom)))
    symmodel = SymbolicModelList(
        Xdom,
        Udom,
        AutomatonList{Set{NTuple{3,Int}}}(0, nu),
        Dict{NTuple{N,T1},Int}(),
        NTuple{N,T1}[],
        upos2int,
        uint2pos,
    )
end


function get_state_by_xpos(
    symmodel::SymbolicModelList{N,M,<:DO.GeneralDomainList{N,DO.RectanglularObstacles{NTuple{N,T}}}},
    pos,
) where {N,M,T}
    id = get(symmodel.xpos2int, pos, nothing)
    created = false
    if id === nothing
        if pos in symmodel.Xdom
            created = true
            push!(symmodel.xint2pos, pos)
            id = length(symmodel.xint2pos)
            symmodel.xpos2int[pos] = id
            i = HybridSystems.add_state!(symmodel.autom)
            @assert i == id
        else
            error("$pos is not in state domain $(symmodel.Xdom)")
        end
    end
    return id::Int,created
end

function get_symbol(symmodel,subset,incl_mode::DO.INCL_MODE)
    Xdom = symmodel.Xdom
    grid = Xdom.grid
    posL = DO.get_subset_pos(Xdom,subset,incl_mode)
    symbolsList = [get_state_by_xpos(symmodel, pos)[1] for pos in posL]
    return symbolsList
end

function get_symbols(symmodel,subsetList,incl_mode::DO.INCL_MODE)
    symbols = Int[]
    for subset in subsetList
        append!(symbols,get_symbol(symmodel,subset,incl_mode))
    end
    return symbols
end
