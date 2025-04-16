"""
    SymbolicModelList{N, M, S1, S2, A} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells.
"""
mutable struct SymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.CustomList{M},
    A,
} <: GridBasedSymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    ucoord2int::Any
    uint2coord::Any
    determinized::Bool
    # Field to store algorithm-specific data
    metadata::Dict{Symbol, Any}
end

function NewSymbolicModelListList(
    Xdom,
    Udom,
    ::Type{S} = UT.SortedTupleSet{3, NTuple{3, Int}},
) where {S}
    nx = DO.get_ncells(Xdom)
    xint2pos = [pos for pos in DO.enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))

    customDomainList = DO.convert_to_custom_domain(Udom)
    nu = DO.get_ncells(customDomainList)
    uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
    ucoord2int =
        Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))
    autom = AutomatonList{S}(nx, nu)

    return SymbolicModelList(
        Xdom,
        customDomainList,
        autom,
        xpos2int,
        xint2pos,
        ucoord2int,
        uint2coord,
        false,
        Dict{Symbol, Any}(),
    )
end

function with_automaton(symmodel::SymbolicModelList, autom)
    return SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        symmodel.ucoord2int,
        symmodel.uint2coord,
        symmodel.determinized,
        symmodel.metadata,
    )
end

get_n_state(symmodel::SymbolicModelList) = length(symmodel.xint2pos)
get_n_input(symmodel::SymbolicModelList) = length(symmodel.uint2coord)
enum_states(symmodel::SymbolicModelList) = 1:get_n_state(symmodel)
enum_inputs(symmodel::SymbolicModelList) = 1:get_n_input(symmodel)
get_state_domain(symmodel::SymbolicModelList) = symmodel.Xdom
get_input_domain(symmodel::SymbolicModelList) = symmodel.Udom

get_xpos_by_state(symmodel::SymbolicModelList, state) = symmodel.xint2pos[state]
get_state_by_xpos(symmodel::SymbolicModelList, xpos) = symmodel.xpos2int[xpos]
is_xpos(symmodel::SymbolicModelList, xpos) = haskey(symmodel.xpos2int, xpos)

function get_transitions(symmodel::SymbolicModelList)
    transitions = get_transition(symmodel.autom)
    return UT.get_data(transitions)
end

function get_concrete_input(symmodel::SymbolicModelList, input)
    if !symmodel.determinized
        return symmodel.uint2coord[input]
    else
        return get_concrete_input(symmodel.metadata[:original_symmodel], symmodel.uint2coord[input][1])
    end
end
function get_abstract_input(symmodel::SymbolicModelList, u)
    if !symmodel.determinized
        return symmodel.ucoord2int[u]
    else
        return get_abstract_input(symmodel.metadata[:original_symmodel], u)
    end
end

add_transitions!(symmodel::SymbolicModelList, translist) =
    add_transitions!(symmodel.autom, translist)

is_deterministic(symmodel::SymbolicModelList) = is_deterministic(symmodel.autom)

function determinize_symbolic_model(symmodel::SymbolicModelList)
    autom = symmodel.autom
    new_uint2coord = Dict{Int, Tuple{Any, Int}}()  # New input mapping (int -> (symbol, target))
    new_ucoord2int = Dict{Tuple{Any, Int}, Int}()  # Reverse (symbol, target) -> int

    next_input_id = 1  # Track new input indices

    # Step 1: Build new input encodings
    for (target, source, symbol) in UT.get_data(get_transition(autom))
        new_input = (symbol, target)  # Couple input with target

        if !haskey(new_ucoord2int, new_input)
            new_ucoord2int[new_input] = next_input_id
            new_uint2coord[next_input_id] = new_input
            next_input_id += 1
        end
    end

     # Step 2: Create new automaton with same number of states, new number of symbols, and transitions
     nstates = autom.nstates
     new_nsymbols = length(new_uint2coord)
     new_autom = AutomatonList{typeof(autom.transitions)}(nstates, new_nsymbols)
     for (target, source, symbol) in UT.get_data(autom.transitions)
        new_input = (symbol, target)
        new_input_id = new_ucoord2int[new_input]
        HybridSystems.add_transition!(new_autom, source, target, new_input_id)
    end

    new_symmodel = SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        new_autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        new_ucoord2int,
        new_uint2coord,
        true,
        Dict{Symbol, Any}(),
    )

    # Step 3: Store the original symmodel to keep the original input mapping (symbol -> u) and (u -> symbol)
    new_symmodel.metadata[:original_symmodel] = symmodel
    return new_symmodel
end
