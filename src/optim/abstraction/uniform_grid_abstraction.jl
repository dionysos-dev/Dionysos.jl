export UniformGridAbstraction

module UniformGridAbstraction

import StaticArrays

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

@enum ApproxMode GROWTH LINEARIZED DELTA_GAS DSICRETE_TIME

using JuMP

"""
    Optimizer{T} <: MOI.AbstractOptimizer

Solver based on the classical abstraction method (used for instance in SCOTS) for which the whole domain is partioned into hyper-rectangular cells, independently of the control task.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, PR.ProblemType}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Any
    discretized_system::Any
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    jacobian_bound::Union{Nothing, Function}
    time_step::T
    num_sub_steps_system_map::Int
    num_sub_steps_growth_bound::Int
    δGAS::Union{Nothing, Bool}
    solve_time_sec::T
    approx_mode::ApproxMode

    ## for the discrete time (no need for the jacobian_bound) 
    #but for the system_map, growthbound_map and sys_inv_map
    system_map::Union{Nothing, Function}
    growthbound_map::Union{Nothing, Function}
    sys_inv_map::Union{Nothing, Function}

    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            NaN,
            5,
            5,
            false,
            0.0,
            GROWTH,
            nothing,
            nothing,
            nothing,
        )
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstraction(concrete_system, model)
    Xfull = DO.DomainList(model.state_grid)
    DO.add_set!(Xfull, concrete_system.X, DO.INNER)
    Ufull = DO.DomainList(model.input_grid)
    DO.add_set!(Ufull, concrete_system.U, DO.CENTER)
    abstract_system = SY.NewSymbolicModelListList(Xfull, Ufull)

    # TODO add noise to the description of the system so in a MathematicalSystems
    #      this is a workaround
    nstates = Dionysos.Utils.get_dims(concrete_system.X)
    noise = StaticArrays.SVector(ntuple(_ -> 0.0, Val(nstates)))

    if isnan(model.time_step)
        error("Please set the `time_step`.")
    end

    model.discretized_system = if model.approx_mode == GROWTH
        if isnothing(model.jacobian_bound)
            error("Please set the `jacobian_bound`.")
        end
        Dionysos.System.NewControlSystemGrowthRK4(
            model.time_step,
            concrete_system.f,
            model.jacobian_bound,
            noise,
            noise,
            model.num_sub_steps_system_map,
            model.num_sub_steps_growth_bound,
        )
    elseif model.approx_mode == DSICRETE_TIME
        if isnothing(model.system_map)
            error("Please set the `system_map`.")
        end
        if isnothing(model.growthbound_map)
            error("Please set the `growthbound_map`.")
        end
        if isnothing(model.sys_inv_map)
            error("Please set the `sys_inv_map`.")
        end
        Dionysos.System.ControlSystemGrowth(
            model.time_step,
            noise,
            noise,
            model.system_map,
            model.growthbound_map,
            model.sys_inv_map,
        )
    elseif model.approx_mode == LINEARIZED
        Dionysos.System.NewControlSystemLinearizedRK4(
            model.time_step,
            concrete_system.f,
            model.jacobian,
            u -> 0.0,
            u -> 0.0,
            noise,
            model.num_sub_steps_system_map,
        )
    else
        @assert model.approx_mode == DELTA_GAS
        Dionysos.System.NewSimpleSystem(
            model.time_step,
            concrete_system.f,
            noise,
            model.num_sub_steps_system_map,
        )
    end

    if model.δGAS
        @time SY.compute_deterministic_symmodel_from_controlsystem!(
            abstract_system,
            model.discretized_system,
        )
    else
        @time SY.compute_symmodel_from_controlsystem!(
            abstract_system,
            model.discretized_system,
        )
    end
    return abstract_system
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xtarget, abstract_system.Xdom, concrete_problem.target_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    target_list =
        [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xtarget)]
    return PR.OptimalControlProblem(
        abstract_system,
        init_list,
        target_list,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xsafe = DO.DomainList(state_grid)
    DO.add_subset!(Xsafe, abstract_system.Xdom, concrete_problem.safe_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    safe_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xsafe)]
    return PR.SafetyProblem(
        abstract_system,
        init_list,
        safe_list,
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function solve_abstract_problem(abstract_problem::PR.OptimalControlProblem)
    abstract_controller = NewControllerList()
    compute_controller_reach!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.target_set,
    )
    return abstract_controller
end

function solve_abstract_problem(abstract_problem::PR.SafetyProblem)
    abstract_controller = NewControllerList()
    compute_controller_safe!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.safe_set,
    )
    return abstract_controller
end

function solve_concrete_problem(abstract_system, abstract_controller)
    function concrete_controller(x; param = false)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        upos = SY.get_upos_by_symbol(abstract_system, symbol)
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u
    end
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Build the abstraction
    abstract_system = build_abstraction(optimizer.concrete_problem.system, optimizer)
    optimizer.abstract_system = abstract_system
    # Build the abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem
    # Solve the abstract problem
    abstract_controller = solve_abstract_problem(abstract_problem)
    optimizer.abstract_controller = abstract_controller
    # Solve the concrete problem
    optimizer.concrete_controller =
        solve_concrete_problem(abstract_system, abstract_controller)

    optimizer.solve_time_sec = time() - t_ref
    return
end

NewControllerList() = UT.SortedTupleSet{2, NTuple{2, Int}}()

function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    for target in 1:(autom.nstates)
        for soursymb in SY.pre(autom, target)
            num_targets_unreachable[soursymb[1], soursymb[2]] += 1
        end
    end
end

function _compute_controller_reach!(
    contr,
    autom,
    init_set,
    target_set,
    num_targets_unreachable,
    current_targets,
    next_targets,
)
    num_init_unreachable = length(init_set)
    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        for target in current_targets
            for (source, symbol) in SY.pre(autom, target)
                if !(source in target_set) &&
                   iszero(num_targets_unreachable[source, symbol] -= 1)
                    push!(target_set, source)
                    push!(next_targets, source)
                    UT.push_new!(contr, (source, symbol))
                    if source in init_set
                        num_init_unreachable -= 1
                    end
                end
            end
        end
        current_targets, next_targets = next_targets, current_targets
    end
    return iszero(num_init_unreachable)
end
function _data(contr, autom, initlist, targetlist)
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    return initset, targetset, num_targets_unreachable, current_targets, next_targets
end
function compute_controller_reach!(contr, autom, initlist, targetlist::Vector{Int})
    println("compute_controller_reach! started")
    # TODO: try to infer whether num_targets_unreachable is sparse or not,
    # and if sparse, use a dictionary instead
    if !_compute_controller_reach!(
        contr,
        autom,
        _data(contr, autom, initlist, targetlist)...,
    )
        println("\ncompute_controller_reach! terminated without covering init set")
        # ProgressMeter.finish!(prog)
        return
    end
    # ProgressMeter.finish!(prog)
    return println("\ncompute_controller_reach! terminated with success")
end

function _compute_pairstable(pairstable, autom)
    for target in 1:(autom.nstates)
        for soursymb in SY.pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end

function compute_controller_safe!(contr, autom, initlist, safelist)
    println("compute_controller_safe! started")
    nstates = autom.nstates
    nsymbols = autom.nsymbols
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]
    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable; dims = 2)
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
        end
    end
    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)
    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end
    nextunsafeset = Set{Int}()

    # prog = ProgressUnknown("# iterations computing controller:")
    while true
        # ProgressMeter.next!(prog)
        for target in unsafeset
            for soursymb in SY.pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end
        if isempty(nextunsafeset)
            break
        end
        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end
    # ProgressMeter.finish!(prog)
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                UT.push_new!(contr, (source, symbol))
            end
        end
    end
    if ⊆(initlist, safeset)
        println("\ncompute_controller_safe! terminated with success")
    else
        println("\ncompute_controller_safe! terminated without covering init set")
    end
end

end
