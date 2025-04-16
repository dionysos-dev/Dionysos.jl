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

function get_transitions(symmodel::SymbolicModel) end
function add_transitions!(symmodel::SymbolicModel, translist) end
function is_deterministic(symmodel::SymbolicModel) end

"""
    GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M}

An intermediate abstract type for symbolic models that rely on a grid-based discretization.
- `N`: Dimension of the state space.
- `M`: Dimension of the input space.
"""
abstract type GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M} end

function get_xpos_by_state(symmodel::GridBasedSymbolicModel, state) end
function get_state_by_xpos(symmodel::GridBasedSymbolicModel, xpos) end
function is_xpos(symmodel::GridBasedSymbolicModel, xpos) end

get_state_grid(symmodel::GridBasedSymbolicModel) = DO.get_grid(get_state_domain(symmodel))

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


@recipe function f(
    symmodel::GridBasedSymbolicModel; 
    arrowsB = false, 
    dims = [1, 2], 
    value_function = [], # Should be a function value_function(state::Int)::Float64 or left as []
    colormap_name = "Blues", 
    default_color = :yellow
)
    # Display the cells
    state_grid = get_state_grid(symmodel)

    if isa(value_function, Function) 
        # 1. Extract needed parts
        projection_map = Dict{Tuple{Int, Int}, Tuple{Float64, Any}}() # value + elem
        for state in enum_states(symmodel)
            v = value_function(state)
            elem = get_concrete_elem(symmodel, state)
            pos = get_xpos_by_state(symmodel, state)
            x1x2 = pos[dims]
            # Only store if it's better (or first time)
            if haskey(projection_map, x1x2)
                if v < projection_map[x1x2][1]
                    projection_map[x1x2] = (v, elem)
                end
            else
                projection_map[x1x2] = (v, elem)
            end
        end
        # 2. Determine maximum finite value (for color scaling)
        finite_vals = filter(isfinite, getindex.(values(projection_map), 1))
        ValueMax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
        # 3. Setup colormap
        cmap = Colors.colormap(colormap_name)
        mycolorMap = UT.Colormap([0.0, ValueMax], cmap)
        # 4. Order states by decreasing value (for proper layering in case of overlapping cells)
        cost_ordered = sort(collect(projection_map); by = x -> -x[2][1])
         # 5. Plot
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
            get_state_domain(symmodel)
        end
    end
    # Display the arrows
    if arrowsB
        for t in get_transitions(symmodel)
            color = RGB(
                abs(0.6 * sin(t[1])),
                abs(0.6 * sin(t[1] + 2π / 3)),
                abs(0.6 * sin(t[1] - 2π / 3)),
            )
            p1 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[2]))
            p2 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[1]))

            @series begin
                color := color
                return t[1] == t[2] ? UT.DrawPoint(p1) : UT.DrawArrow(p1, p2)
            end
        end
    end
end

function compute_abstract_transitions_from_rectangle!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::UT.HyperRectangle,
    abstract_state,
    abstract_input,
    translist,
)
    Xdom = get_state_domain(symmodel)
    rectI = DO.get_pos_lims_outer(DO.get_grid(Xdom), reachable_set)
    ypos_iter = Iterators.product(DO._ranges(rectI)...)
    allin = true
    for ypos in ypos_iter
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = get_state_by_xpos(symmodel, ypos)
        push!(translist, (target, abstract_state, abstract_input))
    end
    return allin
end

function compute_abstract_transitions_from_points!(
    symmodel::GridBasedSymbolicModel,
    reachable_points,
    abstract_state,
    abstract_input,
    translist,
)
    Xdom = get_state_domain(symmodel)
    allin = true
    for y in reachable_points
        ypos = DO.get_pos_by_coord(Xdom, y)
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = get_state_by_xpos(symmodel, ypos)
        push!(translist, (target, abstract_state, abstract_input))
    end
    unique!(translist)
    return allin
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemOverApproximation;
    verbose = false,
    update_interval = Int(1e5),
)
    translist = Tuple{Int, Int, Int}[]
    compute_reachable_set = ST.get_over_approximation_map(concrete_system_approx)
    total_iterations = max(
        div(get_n_input(abstract_system) * get_n_state(abstract_system), update_interval),
        1,
    )
    progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
    count = 0
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        for abstract_state in enum_states(abstract_system)
            concrete_elem = get_concrete_elem(abstract_system, abstract_state)
            reachable_set = compute_reachable_set(concrete_elem, concrete_input)
            empty!(translist)
            allin = compute_abstract_transitions_from_rectangle!(
                abstract_system,
                reachable_set,
                abstract_state,
                abstract_input,
                translist,
            )
            allin && add_transitions!(abstract_system, translist)
            count += 1
            verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
        end
    end
    verbose && ProgressMeter.finish!(progress)
    return
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeGrowthBound;
    verbose = false,
    update_interval = Int(1e5),
)
    translist = Tuple{Int, Int, Int}[]
    growthbound_map = concrete_system_approx.growthbound_map
    system_map = ST.get_system_map(concrete_system_approx)
    r = DO.get_h(DO.get_grid(get_state_domain(abstract_system))) / 2.0
    total_iterations = max(
        div(get_n_input(abstract_system) * get_n_state(abstract_system), update_interval),
        1,
    )
    progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
    count = 0
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fr = growthbound_map(r, concrete_input)
        for abstract_state in enum_states(abstract_system)
            concrete_state = get_concrete_state(abstract_system, abstract_state)
            Fx = system_map(concrete_state, concrete_input)
            reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)
            empty!(translist)
            allin = compute_abstract_transitions_from_rectangle!(
                abstract_system,
                reachable_set,
                abstract_state,
                abstract_input,
                translist,
            )
            allin && add_transitions!(abstract_system, translist)
            count += 1
            verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
        end
    end
    verbose && ProgressMeter.finish!(progress)
    return
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeLinearized;
    verbose = false,
    update_interval = Int(1e5),
)
    Xdom = get_state_domain(abstract_system)
    N = DO.get_dim(Xdom)
    r = DO.get_h(DO.get_grid(Xdom)) / 2.0
    _H_ = SMatrix{N, N}(I) .* r
    _ONE_ = ones(SVector{N})
    e = norm(r, Inf)
    translist = Tuple{Int, Int, Int}[]
    error_map = concrete_system_approx.error_map
    linsys_map = concrete_system_approx.linsys_map
    total_iterations = max(
        div(get_n_input(abstract_system) * get_n_state(abstract_system), update_interval),
        1,
    )
    progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
    count = 0
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fe = error_map(e, concrete_input)
        Fr = r .+ Fe
        for abstract_state in enum_states(abstract_system)
            concrete_state = get_concrete_state(abstract_system, abstract_state)
            Fx, DFx = linsys_map(concrete_state, _H_, concrete_input)
            A = inv(DFx)
            b = abs.(A) * Fr .+ 1.0
            HP = UT.CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = abs.(DFx) * _ONE_ .+ Fe
            reachable_set = UT.HyperRectangle(Fx - rad, Fx + rad)

            empty!(translist)

            allin = compute_abstract_transitions_from_rectangle!(
                abstract_system,
                reachable_set,
                abstract_state,
                abstract_input,
                translist,
            )
            # rectI = DO.get_pos_lims_outer(Xdom.grid, reachable_set)
            # ypos_iter = Iterators.product(DO._ranges(rectI)...)
            # allin = true
            # for ypos in ypos_iter
            #     y = DO.get_coord_by_pos(Xdom.grid, ypos) - Fx
            #     !(y in HP) && continue
            #     if !(ypos in Xdom)
            #         allin = false
            #         break
            #     end
            #     target = get_state_by_xpos(abstract_system, ypos)
            #     push!(translist, (target, concrete_state, abstract_input))
            # end
            allin && add_transitions!(abstract_system, translist)
            count += 1
            verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
        end
    end
    verbose && ProgressMeter.finish!(progress)
    return
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemUnderApproximation;
    verbose = false,
    update_interval = Int(1e5),
)
    translist = Tuple{Int, Int, Int}[]
    under_approximation_map = ST.get_under_approximation_map(concrete_system_approx)
    total_iterations = max(
        div(get_n_input(abstract_system) * get_n_state(abstract_system), update_interval),
        1,
    )
    progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
    count = 0
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        for abstract_state in enum_states(abstract_system)
            concrete_elem = get_concrete_elem(abstract_system, abstract_state)
            reachable_points = under_approximation_map(concrete_elem, concrete_input)
            empty!(translist)
            allin = compute_abstract_transitions_from_points!(
                abstract_system,
                reachable_points,
                abstract_state,
                abstract_input,
                translist,
            )
            allin && add_transitions!(abstract_system, translist)
            count += 1
            verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
        end
    end
    verbose && ProgressMeter.finish!(progress)
    return
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeCenteredSimulation;
    verbose = false,
    update_interval = Int(1e5),
)
    translist = Tuple{Int, Int, Int}[]
    system_map = ST.get_system_map(concrete_system_approx)
    total_iterations = max(
        div(get_n_input(abstract_system) * get_n_state(abstract_system), update_interval),
        1,
    )
    progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
    count = 0
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        for abstract_state in enum_states(abstract_system)
            concrete_state = get_concrete_state(abstract_system, abstract_state)
            reachable_points = [system_map(concrete_state, concrete_input)]
            empty!(translist)
            allin = compute_abstract_transitions_from_points!(
                abstract_system,
                reachable_points,
                abstract_state,
                abstract_input,
                translist,
            )
            allin && add_transitions!(abstract_system, translist)
            count += 1
            verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
        end
    end
    verbose && ProgressMeter.finish!(progress)
    return
end
