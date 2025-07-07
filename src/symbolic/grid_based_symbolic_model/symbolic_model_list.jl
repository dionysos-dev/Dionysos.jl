"""
    SymbolicModelList{N, M, S1, S2, A} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells.
"""
mutable struct SymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.CustomList{M},
    A
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
    AutomatonConstructor::Function = SortedAutomatonListFactory,
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

pre(symmodel::SymbolicModelList, target::Int) = pre(symmodel.autom, target)
post(symmodel::SymbolicModelList, source::Int, input::Int) = post(symmodel.autom, source, input)
enum_transitions(symmodel::SymbolicModelList) = enum_transitions(symmodel.autom)
add_transition!(symmodel::SymbolicModelList, q::Int, q′::Int, u::Int) = add_transition!(symmodel.autom, q, q′, u)
add_transitions!(symmodel::SymbolicModelList, translist) = add_transitions!(symmodel.autom, translist)

function get_concrete_input(symmodel::SymbolicModelList, input)
    if !symmodel.determinized
        return symmodel.uint2coord[input]
    else
        return get_concrete_input(
            symmodel.metadata[:original_symmodel],
            symmodel.uint2coord[input][1],
        )
    end
end
function get_abstract_input(symmodel::SymbolicModelList, u)
    if !symmodel.determinized
        return symmodel.ucoord2int[u]
    else
        return get_abstract_input(symmodel.metadata[:original_symmodel], u)
    end
end


"""
    is_deterministic(symmodel::SymbolicModelList) -> Bool

Returns `true` if the symbolic model is deterministic.

A symbolic model is deterministic if, for every (state, input) pair, there is **at most one successor**.  
"""
is_deterministic(symmodel::SymbolicModelList) = is_deterministic(symmodel.autom)

"""
    determinize_symbolic_model(symmodel::SymbolicModelList) -> SymbolicModelList

Returns a **determinized version** of the given symbolic model by encoding each transition input as a pair `(input_symbol, target_state)`.

This transformation removes nondeterminism by lifting the input space: each original input is paired with its intended target, making transitions unique.  

# Arguments
- `symmodel`: A [`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList) that may contain nondeterministic transitions.

# Returns
- A new [`SymbolicModelList`](@ref) that is deterministic, with:
    - New unique inputs: `(symbol, target)` pairs.
    - Transitions remapped to be deterministic.
    - `metadata[:original_symmodel]` containing a reference to the original model.

# Notes
- The determinized model will have more inputs than the original.

# Example

```julia
det_symmodel = determinize_symbolic_model(symmodel)
is_deterministic(det_symmodel) == true
```
"""
function determinize_symbolic_model(symmodel::SymbolicModelList; AutomatonConstructor::Function = SortedAutomatonListFactory)
    transitions = enum_transitions(symmodel)

    new_ucoord2int = Dict{Tuple{Any, Int}, Int}() # New input mapping (int -> (symbol, target))
    new_uint2coord = Dict{Int, Tuple{Any, Int}}() # Reverse (symbol, target) -> int

    input_counter = 1
    transition_buffer = Vector{NTuple{3, Int}}()
    for (target, source, symbol) in transitions
        new_input = (symbol, target)

        # Get or assign a new input ID
        input_id = get!(new_ucoord2int, new_input) do
            new_uint2coord[input_counter] = new_input
            input_counter += 1
            return input_counter - 1
        end

        push!(transition_buffer, (target, source, input_id))
    end
    new_autom = AutomatonConstructor(get_n_state(symmodel), length(new_uint2coord))
    add_transitions!(new_autom, transition_buffer)

    # Build new symbolic model
    new_symmodel = SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        new_autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        new_ucoord2int,
        new_uint2coord,
        true,  # mark as determinized
        Dict{Symbol, Any}(),
    )

    new_symmodel.metadata[:original_symmodel] = symmodel

    return new_symmodel
end
