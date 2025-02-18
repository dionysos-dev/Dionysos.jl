# Abstract system where each state contain itself an abstract system.
mutable struct HierarchicalSymbolicSystem # <: SymbolicModel
    symmodel::Any
    sub_symmodels::Any # Dictionnary: (Int, symmodel)
end

function get_sub_symmodel(hierarchical_symmodel::HierarchicalSymbolicSystem, state)
    return hierarchical_symmodel.sub_symmodels[state]
end

function get_xpos_by_state(hierarchical_symmodel::HierarchicalSymbolicSystem, state)
    return get_xpos_by_state(hierarchical_symmodel.symmodel, state)
end

function get_state_by_xpos(hierarchical_symmodel::HierarchicalSymbolicSystem, pos)
    return get_state_by_xpos(hierarchical_symmodel.symmodel, pos)
end

function get_state_by_coord(hierarchical_symmodel::HierarchicalSymbolicSystem, coord)
    return get_state_by_coord(hierarchical_symmodel.symmodel, coord)
end

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
    return get_states_from_sets(hierarchical_symmodel.symmodel, subsetList, incl_mode)
end

function get_ncells(hierarchical_symmodel::HierarchicalSymbolicSystem)
    return get_ncells(hierarchical_symmodel.symmodel)
end

function enum_states(hierarchical_symmodel::HierarchicalSymbolicSystem)
    return enum_states(hierarchical_symmodel.symmodel)
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
