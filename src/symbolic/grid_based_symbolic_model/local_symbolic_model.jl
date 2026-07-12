# ------------------------------------------------------------
# Local partition symbolic model
# ------------------------------------------------------------

"""
    LocalGridBasedSymbolicModel

Wrapper around a global symbolic model that overrides only the source domain.
The state mapping, input mapping, retained domain and input domain remain global.
"""
struct LocalGridBasedSymbolicModel{N, M, SM, XS} <: GridBasedSymbolicModel{N, M}
    parent::SM
    Xset_local::XS
end

function LocalGridBasedSymbolicModel(
    symmodel::GridBasedSymbolicModel{N, M},
    Xset_local,
) where {N, M}
    return LocalGridBasedSymbolicModel{N, M, typeof(symmodel), typeof(Xset_local)}(
        symmodel,
        Xset_local,
    )
end

get_state_mapping(sym::LocalGridBasedSymbolicModel) = get_state_mapping(sym.parent)

get_input_mapping(sym::LocalGridBasedSymbolicModel) = get_input_mapping(sym.parent)

get_state_set(sym::LocalGridBasedSymbolicModel) = sym.Xset_local

get_retained_set(sym::LocalGridBasedSymbolicModel) = get_retained_set(sym.parent)

get_input_set(sym::LocalGridBasedSymbolicModel) = get_input_set(sym.parent)

get_concrete_state(sym::LocalGridBasedSymbolicModel, q) = get_concrete_state(sym.parent, q)

get_concrete_elem(sym::LocalGridBasedSymbolicModel, q::Int) =
    get_concrete_elem(sym.parent, q)

get_concrete_input(sym::LocalGridBasedSymbolicModel, u) = get_concrete_input(sym.parent, u)

get_abstract_state(sym::LocalGridBasedSymbolicModel, x) = get_abstract_state(sym.parent, x)

add_transitions!(sym::LocalGridBasedSymbolicModel, trans) =
    add_transitions!(sym.parent, trans)

get_transition_metadata(sym::LocalGridBasedSymbolicModel) =
    get_transition_metadata(sym.parent)
