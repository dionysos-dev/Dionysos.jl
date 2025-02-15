using Plots, Colors

"""
    Abstract Type: SymbolicModel{N, M}

Defines a generic symbolic model interface, where:
- `N` is the state space dimension.
- `M` is the input space dimension.
"""
abstract type SymbolicModel{N, M} end

function get_n_state(symmodel::SymbolicModel) end
function get_n_input(symmodel::SymbolicModel) end
function enum_states(symmodel::SymbolicModel) end
function enum_inputs(symmodel::SymbolicModel) end
function get_state_domain(symmodel::SymbolicModel) end
function get_input_domain(symmodel::SymbolicModel) end

function get_concrete_state(symmodel::SymbolicModel, state) end
function get_concrete_input(symmodel::SymbolicModel, input) end
function get_abstract_state(symmodel::SymbolicModel, x) end
function get_abstract_input(symmodel::SymbolicModel, u) end

function add_transitions!(symmodel::SymbolicModel, translist) end

"""
    GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M}

An intermediate abstract type for symbolic models that rely on a grid-based discretization.
- `N`: Dimension of the state space.
- `M`: Dimension of the input space.
"""
abstract type GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M} end

function get_xpos_by_state(symmodel::GridBasedSymbolicModel, state) end
function get_state_by_xpos(symmodel::GridBasedSymbolicModel, xpos) end
# function get_upos_by_symbol(symmodel::GridBasedSymbolicModel, symbol) end
# function get_symbol_by_upos(symmodel::GridBasedSymbolicModel, upos) end
function is_xpos(symmodel::GridBasedSymbolicModel, xpos) end

function get_concrete_state(symmodel::GridBasedSymbolicModel, state)
    xpos = get_xpos_by_state(symmodel, state)
    return DO.get_coord_by_pos(get_state_domain(symmodel), xpos)
end

function get_concrete_elem(symmodel::GridBasedSymbolicModel, state)
    xpos = get_xpos_by_state(symmodel, state)
    return DO.get_elem_by_pos(get_state_domain(symmodel), xpos)
end

function get_abstract_state(symmodel::GridBasedSymbolicModel, x)
    xpos = DO.get_pos_by_coord(get_state_domain(symmodel), x)
    return get_state_by_xpos(symmodel, xpos)
end

get_state_grid(symmodel::GridBasedSymbolicModel) = DO.get_grid(get_state_domain(symmodel))
get_states_by_xpos(symmodel::GridBasedSymbolicModel, l_xpos) =
    [get_state_by_xpos(symmodel, xpos) for xpos in l_xpos]

function get_domain_from_states(symmodel::GridBasedSymbolicModel, states)
    newDomain = DO.DomainList(get_state_grid(symmodel))
    for state in states
        DO.add_pos!(newDomain, get_xpos_by_state(symmodel, state))
    end
    return newDomain
end

function get_states_from_set(
    symmodel::GridBasedSymbolicModel,
    subset::UT.HyperRectangle,
    incl_mode::DO.INCL_MODE,
)
    Xdom = get_state_domain(symmodel)
    posL = DO.get_subset_pos(Xdom, subset, incl_mode)
    return [get_state_by_xpos(symmodel, pos) for pos in posL]
end

function get_states_from_sets(
    symmodel::GridBasedSymbolicModel,
    subsets,
    incl_mode::DO.INCL_MODE,
)
    states = []
    for subset in subsets
        append!(states, get_states_from_set(symmodel, subset, incl_mode))
    end
    return states
end
