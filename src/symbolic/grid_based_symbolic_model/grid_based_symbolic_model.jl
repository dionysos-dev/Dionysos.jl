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

using ProgressMeter
function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeCenteredSimulation;
    verbose = false,
    update_interval = Int(1e3),
)
    translist = Tuple{Int, Int, Int}[]
    system_map = ST.get_system_map(concrete_system_approx)
    n_forward_images = get_n_input(abstract_system) * get_n_state(abstract_system)
    total_iterations = max(
        div(n_forward_images, update_interval),
        1,
    )
    iteration_progress_percentage = update_interval / n_forward_images
    @info("Total number of bar updates : $total_iterations, the progress bar will update every $(iteration_progress_percentage)%")
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
