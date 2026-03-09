"""
SymbolicModelList:
- XMapping / UMapping define universes and conversions
- Xset: states we build from (sources)
- Rset: allowed states as targets (superset; "relation universe allowance")
- Uset: inputs considered
"""
mutable struct SymbolicModelList{
    N,
    M,
    TX,
    TU,
    XM <: MP.AbstractMapping{N, TX},
    UM <: MP.AbstractMapping{M, TU},
    XS <: MP.AbstractStateSet{N},
    RS <: MP.AbstractStateSet{N},
    US <: MP.AbstractStateSet{M},
    A,
    OS,
} <: GridBasedSymbolicModel{N, M}
    XMapping::XM
    UMapping::UM
    Xset::XS
    Rset::RS
    Uset::US
    autom::A
    original_symmodel::OS
end

# default sets = "all states of mapping"
_default_stateset(::MP.AbstractMapping{N, TX}) where {N, TX} = MP.MappingSet{N}()

function SymbolicModelList(
    XMapping::XM,
    UMapping::UM;
    Xset::Union{Nothing, MP.AbstractStateSet{N}} = nothing,
    Rset::Union{Nothing, MP.AbstractStateSet{N}} = nothing,
    Uset::Union{Nothing, MP.AbstractStateSet{M}} = nothing,
    automaton_constructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
    original_symmodel = nothing,
    convert_U_to_list::Bool = true,
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
        XMapping,
        UMap,
        Xset_final,
        Rset_final,
        Uset_final,
        autom,
        original_symmodel,
    )
end

# --- interface ---
get_state_mapping(sym::SymbolicModelList) = sym.XMapping
get_input_mapping(sym::SymbolicModelList) = sym.UMapping

# get_state_domain(sym::SymbolicModelList) = sym.Xset
get_source_domain(sym::SymbolicModelList) = sym.Xset
# get_retained_domain(sym::SymbolicModelList) = sym.Rset
get_retained_domain(sym::SymbolicModelList) = sym.Rset
get_input_domain(sym::SymbolicModelList) = sym.Uset

get_automaton(sym::SymbolicModelList) = sym.autom

is_determinized(sym::SymbolicModelList) = !(sym.original_symmodel === nothing)
