"""
    SymbolicModelList{N, M, S1, S2, A, U, OS} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells.
"""
mutable struct SymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.CustomList{M},
    A,
    U,
    OS,
} <: GridBasedSymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    ucoord2int::Dict{U, Int}
    uint2coord::Vector{U}
    original_symmodel::OS
end

"""
    SymbolicModelList(Xdom, Udom; AutomatonConstructor)

Constructor for a fresh (non-determinized) SymbolicModelList.
"""
function SymbolicModelList(
    Xdom,
    Udom,
    AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
)
    nx = DO.get_ncells(Xdom)
    xint2pos = [pos for pos in DO.enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))

    customDomainList = DO.convert_to_custom_domain(Udom)
    nu = DO.get_ncells(customDomainList)
    uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
    ucoord2int =
        Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))
    autom = AutomatonConstructor(nx, nu)

    return SymbolicModelList(
        Xdom,
        customDomainList,
        autom,
        xpos2int,
        xint2pos,
        ucoord2int,
        uint2coord,
        nothing, # OS = Nothing for base models
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
        symmodel.original_symmodel,
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

pre(symmodel::SymbolicModelList, target::Int) = pre(symmodel.autom, target)
post(symmodel::SymbolicModelList, source::Int, input::Int) =
    post(symmodel.autom, source, input)
enum_transitions(symmodel::SymbolicModelList) = enum_transitions(symmodel.autom)
add_transition!(symmodel::SymbolicModelList, q::Int, q′::Int, u::Int) =
    add_transition!(symmodel.autom, q, q′, u)
add_transitions!(symmodel::SymbolicModelList, translist) =
    add_transitions!(symmodel.autom, translist)

is_determinized(symmodel::SymbolicModelList) = !(symmodel.original_symmodel === nothing)

"""
    is_deterministic(symmodel::SymbolicModelList) -> Bool

Returns `true` if the symbolic model is deterministic.
"""
is_deterministic(symmodel::SymbolicModelList) = is_deterministic(symmodel.autom)

function get_concrete_input(
    symmodel::SymbolicModelList{N, M, S1, S2, A, U, Nothing},
    input::Int,
) where {N, M, S1, S2, A, U}
    return symmodel.uint2coord[input]
end

function get_abstract_input(
    symmodel::SymbolicModelList{N, M, S1, S2, A, U, Nothing},
    u::U,
) where {N, M, S1, S2, A, U}
    return symmodel.ucoord2int[u]
end

function get_concrete_input(
    symmodel::SymbolicModelList{N, M, S1, S2, A, Tuple{Uprev, Int}, OS},
    input::Int,
) where {N, M, S1, S2, A, Uprev, OS}
    u, _ = symmodel.uint2coord[input]
    return get_concrete_input(symmodel.original_symmodel, u)
end

function get_abstract_input(
    symmodel::SymbolicModelList{N, M, S1, S2, A, Tuple{Uprev, Int}, OS},
    u,
) where {N, M, S1, S2, A, Uprev, OS}
    return get_abstract_input(symmodel.original_symmodel, u)
end

"""
    determinize_symbolic_model(symmodel::SymbolicModelList) -> SymbolicModelList

Returns a determinized version of the given symbolic model by encoding each transition input as a pair `(input_symbol, target_state)`.
"""
function determinize_symbolic_model(
    symmodel::SymbolicModelList;
    AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
)
    transitions = enum_transitions(symmodel)

    U = eltype(symmodel.uint2coord)
    new_ucoord2int = Dict{Tuple{U, Int}, Int}()
    new_uint2coord = Tuple{U, Int}[]

    transition_buffer = Vector{NTuple{3, Int}}()
    for (target, source, symbol) in transitions
        u_coord = symmodel.uint2coord[symbol]  # Get symbolic input
        new_input = (u_coord, target)          # Determinize with symbolic input and target state
        input_id = get!(new_ucoord2int, new_input) do
            push!(new_uint2coord, new_input)
            return length(new_uint2coord)
        end

        push!(transition_buffer, (target, source, input_id))
    end
    new_autom = AutomatonConstructor(get_n_state(symmodel), length(new_uint2coord))
    add_transitions!(new_autom, transition_buffer)

    new_symmodel = SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        new_autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        new_ucoord2int,
        new_uint2coord,
        symmodel,
    )

    return new_symmodel
end
