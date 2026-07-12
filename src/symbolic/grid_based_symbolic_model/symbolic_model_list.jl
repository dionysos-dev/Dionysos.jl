"""
SymbolicModelList:
- `X`: source states bundled with the state mapping (`MP.MappedStateSet`)
- `R`: allowed target states, same state mapping ("relation universe allowance")
- `U`: inputs considered, bundled with the input mapping
"""
mutable struct SymbolicModelList{
    N,
    M,
    XB <: MP.MappedStateSet{N},
    RB <: MP.MappedStateSet{N},
    UB <: MP.MappedStateSet{M},
    A,
    OS,
    MD <: AbstractTransitionMetadata,
} <: GridBasedSymbolicModel{N, M}
    X::XB
    R::RB
    U::UB
    autom::A
    original_symmodel::OS
    metadata::MD
end

# default sets = "all states of mapping"
_default_stateset(::MP.AbstractMapping{N, TX}) where {N, TX} = MP.FullStateSet{N}()

function SymbolicModelList(
    XMapping::XM,
    UMapping::UM;
    Xset::Union{Nothing, MP.AbstractStateSet{N}} = nothing,
    Rset::Union{Nothing, MP.AbstractStateSet{N}} = nothing,
    Uset::Union{Nothing, MP.AbstractStateSet{M}} = nothing,
    automaton_constructor::Function = (n, m) -> SortedAutomatonList(n, m),
    original_symmodel = nothing,
    convert_U_to_list::Bool = true,
    metadata = NoTransitionMetadata(),
) where {N, M, TX, TU, XM <: MP.AbstractMapping{N, TX}, UM <: MP.AbstractMapping{M, TU}}
    UMap = convert_U_to_list ? MP.convert_to_list_mapping(UMapping) : UMapping

    Xset_final = Xset === nothing ? _default_stateset(XMapping) : Xset
    Uset_final = Uset === nothing ? _default_stateset(UMap) : Uset
    Rset_final = Rset === nothing ? Xset_final : Rset

    autom = automaton_constructor(
        maximum(MP.enum_states(Rset_final, XMapping)), #MP.get_n_state(Xset_final, XMapping),
        MP.get_n_state(Uset_final, UMap),
    )

    return SymbolicModelList(
        MP.MappedStateSet(Xset_final, XMapping),
        MP.MappedStateSet(Rset_final, XMapping),
        MP.MappedStateSet(Uset_final, UMap),
        autom,
        original_symmodel,
        metadata,
    )
end

# --- interface ---
get_mapped_state_set(sym::SymbolicModelList) = sym.X
get_mapped_retained_set(sym::SymbolicModelList) = sym.R
get_mapped_input_set(sym::SymbolicModelList) = sym.U

get_state_mapping(sym::SymbolicModelList) = MP.get_mapping(sym.X)
get_input_mapping(sym::SymbolicModelList) = MP.get_mapping(sym.U)

get_state_set(sym::SymbolicModelList) = MP.get_set(sym.X)
get_retained_set(sym::SymbolicModelList) = MP.get_set(sym.R)
get_input_set(sym::SymbolicModelList) = MP.get_set(sym.U)

get_automaton(sym::SymbolicModelList) = sym.autom

is_determinized(sym::SymbolicModelList) = !(sym.original_symmodel === nothing)

get_transition_metadata(sym::SymbolicModelList) = sym.metadata

function without_metadata(sym::SymbolicModel)
    return SymbolicModelList(
        get_state_mapping(sym),
        get_input_mapping(sym);
        Xset = get_state_set(sym),
        Rset = get_retained_set(sym),
        Uset = get_input_set(sym),
        automaton_constructor = (n, m) -> get_automaton(sym),
        metadata = NoTransitionMetadata(),
    )
end
