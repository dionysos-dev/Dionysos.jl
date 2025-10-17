abstract type AbstractRelation end

get_concrete_system(::AbstractRelation)
get_abstract_system(::AbstractRelation)

# Mappings between concrete and abstract states
get_abstract_states(::AbstractRelation, x::ConcreteState)
get_concrete_states(::AbstractRelation, s::Int)

# Convenience for partition-based relations
get_abstract_state(::AbstractRelation, x::ConcreteState)
get_concrete_state(::AbstractRelation, s::Int)

# Input mappings
get_abstract_input(::AbstractRelation, u::ConcreteInput)
get_concrete_input(::AbstractRelation, symb::Int)

# Domain info
get_state_domain(::AbstractRelation)
get_input_domain(::AbstractRelation)

# Boolean utilities
is_deterministic(::AbstractRelation)
is_abstract_system_deterministic(R::AbstractRelation) = is_deterministic(get_abstract_system(R))

# Size/structure queries
get_n_state_abstract_system(R::AbstractRelation) = SY.get_n_state(get_abstract_system(R))
get_n_input_abstract_system(R::AbstractRelation) = SY.get_n_input(get_abstract_system(R))

# Return the abstraction type :FRR, :ASR, :MCR, etc.
abstraction_type(::AbstractRelation)



############################################
################ Grid Based ################
############################################

abstract type GridBasedRelation <: AbstractRelation end

get_state_grid(R::GridBasedRelation) = DO.get_grid(get_state_domain(R))

get_state_by_pos(::GridBasedRelation, pos::NTuple)
get_pos_by_state(::GridBasedRelation, s::Int)
is_pos(symmodel::GridBasedRelation, pos)

get_state_grid(symmodel::GridBasedRelation) = DO.get_grid(get_state_domain(symmodel))
function get_concrete_state(symmodel::GridBasedRelation, state)
    pos = get_state_by_pos(symmodel, state)
    return DO.get_coord_by_pos(get_state_domain(symmodel), pos)
end

function get_concrete_elem(symmodel::GridBasedRelation, state)
    pos = get_state_by_pos(symmodel, state)
    return DO.get_elem_by_pos(get_state_domain(symmodel), pos)
end

function get_abstract_state(symmodel::GridBasedRelation, x)
    pos = DO.get_pos_by_coord(get_state_domain(symmodel), x)
    return get_state_by_pos(symmodel, pos)
end

get_states_by_pos(symmodel::GridBasedRelation, l_pos) =
    [get_state_by_pos(symmodel, pos) for pos in l_pos]

function get_domain_from_states(symmodel::GridBasedRelation, states)
    newDomain = DO.DomainList(get_state_grid(symmodel))
    for state in states
        DO.add_pos!(newDomain, get_state_by_pos(symmodel, state))
    end
    return newDomain
end

function get_states_from_set(
    symmodel::GridBasedRelation,
    subset::UT.HyperRectangle,
    incl_mode::DO.INCL_MODE,
)
    Xdom = get_state_domain(symmodel)
    posL = DO.get_subset_pos(Xdom, subset, incl_mode)
    return [get_state_by_pos(symmodel, pos) for pos in posL]
end

function get_states_from_sets(
    symmodel::GridBasedRelation,
    subsets::UT.LazyUnionSetArray,
    incl_mode::DO.INCL_MODE,
)
    states = []
    for subset in subsets.sets
        append!(states, get_states_from_set(symmodel, subset, incl_mode))
    end
    return states
end

function get_states_from_sets(
    symmodel::GridBasedRelation,
    subsets,
    incl_mode::DO.INCL_MODE,
)
    states = []
    for subset in subsets
        append!(states, get_states_from_set(symmodel, subset, incl_mode))
    end
    return states
end

@recipe function f(
    r::GridBasedRelation;
    arrowsB = false,
    dims = [1, 2],
    value_function = [],
    colormap_name = "Blues",
    default_color = :yellow,
)
    state_grid = get_state_grid(r)
    projection_map = Dict{Tuple{Int, Int}, Tuple{Float64, Any}}()

    if isa(value_function, Function)
        for (q, pos) in enumerate(r.state2pos)
            v = value_function(q)
            elem = DO.get_elem_by_pos(r.Xdom, pos)
            x1x2 = pos[dims]
            if haskey(projection_map, x1x2)
                if v < projection_map[x1x2][1]
                    projection_map[x1x2] = (v, elem)
                end
            else
                projection_map[x1x2] = (v, elem)
            end
        end

        finite_vals = filter(isfinite, getindex.(values(projection_map), 1))
        ValueMax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
        cmap = Colors.colormap(colormap_name)
        mycolorMap = UT.Colormap([0.0, ValueMax], cmap)
        cost_ordered = sort(collect(projection_map); by = x -> -x[2][1])

        @series begin
            for (_, (value, elem)) in cost_ordered
                color = isfinite(value) ? UT.get_color(mycolorMap, value) : default_color
                @series begin
                    color := color
                    dims := dims
                    label := ""
                    return elem
                end
            end
            mycolorMap
        end
    else
        @series begin
            dims := dims
            r.Xdom
        end
    end
end




