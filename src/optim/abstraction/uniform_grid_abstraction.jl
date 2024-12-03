export UniformGridAbstraction

module UniformGridAbstraction

import StaticArrays

import MathematicalSystems

import Dionysos

@enum ApproxMode GROWTH LINEARIZED DELTA_GAS

using JuMP

"""
    Optimizer{T} <: MOI.AbstractOptimizer

Solver based on the classical abstraction method (used for instance in SCOTS) for which the whole domain is partioned into hyper-rectangular cells, independently of the control task.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, Dionysos.Problem.ProblemType}
    abstract_problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem, Dionysos.Problem.SafetyProblem}
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}
    abstract_controller::Union{Nothing, Dionysos.Utils.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Any
    discretized_system::Any
    state_grid::Union{Nothing, Dionysos.Domain.Grid}
    input_grid::Union{Nothing, Dionysos.Domain.Grid}
    jacobian_bound::Union{Nothing, Function}
    time_step::T
    num_sub_steps_system_map::Int
    num_sub_steps_growth_bound::Int
    δGAS::Union{Nothing, Bool}
    solve_time_sec::T
    approx_mode::ApproxMode

    ## for the discrete time (no need for the jacobian_bound) 
    # but for the `growthbound_map` and `sys_inv_map`
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

function discretized_system(
    concrete_system::MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem,
    model,
    noise,
)
    if isnothing(model.growthbound_map)
        error("Please set the `growthbound_map`.")
    end
    if isnothing(model.sys_inv_map)
        error("Please set the `sys_inv_map`.")
    end
    return Dionysos.System.ControlSystemGrowth(
        1.0, # `time_step` should be ignored by `concrete_system.f`, `model.growthbound_map` and `model.sys_inv_map` anyway
        noise,
        noise,
        concrete_system.f,
        model.growthbound_map,
        model.sys_inv_map,
    )
end

function discretized_system(
    concrete_system::MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem,
    model,
    noise,
)
    if isnan(model.time_step)
        error("Please set the `time_step`.")
    end

    if model.approx_mode == GROWTH
        if isnothing(model.jacobian_bound)
            error("Please set the `jacobian_bound`.")
        end
        return Dionysos.System.discretize_system_with_growth_bound(
            model.time_step,
            concrete_system.f,
            model.jacobian_bound,
            noise,
            noise,
            model.num_sub_steps_system_map,
            model.num_sub_steps_growth_bound,
        )
    elseif model.approx_mode == LINEARIZED
        return Dionysos.System.discretize_system_with_linearization(
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
        return Dionysos.System.NewSimpleSystem(
            model.time_step,
            concrete_system.f,
            noise,
            model.num_sub_steps_system_map,
        )
    end
end

function build_abstraction(concrete_system, model)
    if isnothing(model.state_grid)
        error("Please set the `state_grid`.")
    end

    if isnothing(model.input_grid)
        error("Please set the `input_grid`.")
    end

    Xfull = Dionysos.Domain.DomainList(model.state_grid)
    Dionysos.Domain.add_set!(Xfull, concrete_system.X, Dionysos.Domain.INNER)
    Ufull = Dionysos.Domain.DomainList(model.input_grid)
    Dionysos.Domain.add_set!(Ufull, concrete_system.U, Dionysos.Domain.CENTER)
    abstract_system = Dionysos.Symbolic.NewSymbolicModelListList(Xfull, Ufull)

    # TODO add noise to the description of the system so in a MathematicalSystems
    #      this is a workaround
    nstates = Dionysos.Utils.get_dims(concrete_system.X)
    noise = StaticArrays.SVector(ntuple(_ -> 0.0, Val(nstates)))

    model.discretized_system = discretized_system(concrete_system, model, noise)

    if model.δGAS
        @time Dionysos.Symbolic.compute_deterministic_symmodel_from_controlsystem!(
            abstract_system,
            model.discretized_system,
        )
    else
        @time Dionysos.Symbolic.compute_symmodel_from_controlsystem!(
            abstract_system,
            model.discretized_system,
        )
    end
    return abstract_system
end

function build_abstract_problem(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, Dionysos.Domain.OUTER)
    Xtarget = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xtarget, abstract_system.Xdom, concrete_problem.target_set, Dionysos.Domain.INNER)
    init_list = [Dionysos.Symbolic.get_state_by_xpos(abstract_system, pos) for pos in Dionysos.Domain.enum_pos(Xinit)]
    target_list =
        [Dionysos.Symbolic.get_state_by_xpos(abstract_system, pos) for pos in Dionysos.Domain.enum_pos(Xtarget)]
    return Dionysos.Problem.OptimalControlProblem(
        abstract_system,
        init_list,
        target_list,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    concrete_problem::Dionysos.Problem.SafetyProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, Dionysos.Domain.OUTER)
    Xsafe = Dionysos.Domain.DomainList(state_grid)
    Dionysos.Domain.add_subset!(Xsafe, abstract_system.Xdom, concrete_problem.safe_set, Dionysos.Domain.INNER)
    init_list = [Dionysos.Symbolic.get_state_by_xpos(abstract_system, pos) for pos in Dionysos.Domain.enum_pos(Xinit)]
    safe_list = [Dionysos.Symbolic.get_state_by_xpos(abstract_system, pos) for pos in Dionysos.Domain.enum_pos(Xsafe)]
    return Dionysos.Problem.SafetyProblem(
        abstract_system,
        init_list,
        safe_list,
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function solve_abstract_problem(abstract_problem::Dionysos.Problem.OptimalControlProblem)
    abstract_controller = NewControllerList()
    compute_controller_reach!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.target_set,
    )
    return abstract_controller
end

function solve_abstract_problem(abstract_problem::Dionysos.Problem.SafetyProblem)
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
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
        symbollist = Dionysos.Utils.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        upos = Dionysos.Symbolic.get_upos_by_symbol(abstract_system, symbol)
        u = Dionysos.Domain.get_coord_by_pos(abstract_system.Udom.grid, upos)
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

NewControllerList() = Dionysos.Utils.SortedTupleSet{2, NTuple{2, Int}}()

function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    for target in 1:(autom.nstates)
        for soursymb in Dionysos.Symbolic.pre(autom, target)
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
            for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
                if !(source in target_set) &&
                   iszero(num_targets_unreachable[source, symbol] -= 1)
                    push!(target_set, source)
                    push!(next_targets, source)
                    Dionysos.Utils.push_new!(contr, (source, symbol))
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
        for soursymb in Dionysos.Symbolic.pre(autom, target)
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
            for soursymb in Dionysos.Symbolic.pre(autom, target)
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
                Dionysos.Utils.push_new!(contr, (source, symbol))
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
