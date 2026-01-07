# Abstract system where each state contain itself an abstract system.
mutable struct HierarchicalSymbolicSystem # <: SymbolicModel
    symmodel::Any
    sub_symmodels::Any # Dictionnary: (Int, symmodel)
end

get_n_state(hierarchical_symmodel::HierarchicalSymbolicSystem) =
    get_n_state(hierarchical_symmodel.symmodel)
get_n_input(hierarchical_symmodel::HierarchicalSymbolicSystem) =
    get_n_input(hierarchical_symmodel.symmodel)
enum_states(hierarchical_symmodel::HierarchicalSymbolicSystem) =
    enum_states(hierarchical_symmodel.symmodel)
enum_inputs(hierarchical_symmodel::HierarchicalSymbolicSystem) =
    enum_inputs(hierarchical_symmodel.symmodel)

get_xpos_by_state(hierarchical_symmodel::HierarchicalSymbolicSystem, state) =
    get_xpos_by_state(hierarchical_symmodel.symmodel, state)
get_state_by_xpos(hierarchical_symmodel::HierarchicalSymbolicSystem, pos) =
    get_state_by_xpos(hierarchical_symmodel.symmodel, pos)
get_state_by_coord(hierarchical_symmodel::HierarchicalSymbolicSystem, coord) =
    get_state_by_coord(hierarchical_symmodel.symmodel, coord)

get_sub_symmodel(hierarchical_symmodel::HierarchicalSymbolicSystem, state) =
    hierarchical_symmodel.sub_symmodels[state]

function get_symbol(
    hierarchical_symmodel::HierarchicalSymbolicSystem,
    subset,
    incl_mode::DO.INCL_MODE,
)
    return get_states_from_set(hierarchical_symmodel.symmodel, subset, incl_mode)
end

function get_symbols(
    hierarchical_symmodel::HierarchicalSymbolicSystem,
    subsetList,
    incl_mode::DO.INCL_MODE,
)
    return get_states_from_set(hierarchical_symmodel.symmodel, subsetList, incl_mode)
end

@recipe function f(
    hierarchical_symmodel::HierarchicalSymbolicSystem;
    dims = [1, 2],
    arrowsB = true,
    cost = false,
    lyap_fun = [],
)
    @series begin
        arrowsB := false
        hierarchical_symmodel.symmodel
    end
    for (state, sub_symmodel) in hierarchical_symmodel.sub_symmodels
        @series begin
            sub_symmodel
        end
    end
end
