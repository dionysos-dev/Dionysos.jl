# ============================================================
# Optimal Control
# ============================================================

mutable struct OptimizerOptimalControlProblem{T} <: MOI.AbstractOptimizer
    # inputs
    problem::Union{Nothing, PR.OptimalControlProblem}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # outputs
    controller::Union{Nothing, MS.ConstrainedBlackBoxMap}
    controllable_set::Any
    uncontrollable_set::Any
    value_fun_tab::Any
    value_function::Any
    success::Bool
    solve_time_sec::T

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing, # problem
            false,   # early_stop
            false,   # sparse_input
            1,       # print_level
            nothing, # controller
            nothing, # controllable_set
            nothing, # uncontrollable_set
            nothing, # value_fun_tab
            nothing, # value_function
            false,   # success
            zero(T), # solve_time_sec
        )
    end
end

OptimizerOptimalControlProblem() = OptimizerOptimalControlProblem{Float64}()

MOI.is_empty(optimizer::OptimizerOptimalControlProblem) = optimizer.problem === nothing

function MOI.set(
    model::OptimizerOptimalControlProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerOptimalControlProblem, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function MOI.get(model::OptimizerOptimalControlProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_value_function(value_fun_tab)
    return state -> value_fun_tab[state]
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system
    init_set = optimizer.early_stop ? problem.initial_set : ST.enum_states(autom)

    controller, controllable_set, uncontrollable_set, value_fun_tab =
        compute_worst_case_cost_controller(
            autom,
            problem.target_set;
            initial_set = init_set,
            sparse_input = optimizer.sparse_input,
            cost_function = problem.transition_cost,
        )

    optimizer.controller = controller
    optimizer.controllable_set = controllable_set
    optimizer.uncontrollable_set = uncontrollable_set
    optimizer.value_fun_tab = value_fun_tab
    optimizer.value_function = build_value_function(value_fun_tab)
    optimizer.success = all(q -> q in controllable_set, problem.initial_set)

    optimizer.print_level >= 1 &&
        println("Optimal control terminated with success = $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end

"""
    compute_worst_case_cost_controller(
        autom::ST.AbstractAutomatonList,
        target_set;
        initial_set = ST.enum_states(autom),
        cost_function = nothing,
        sparse_input::Bool = false,
    )

Compute controller for optimal control problem assuming **worst case** for uncertainty.
This means, if there are several states `x⁺` that can be the next state from state
`x` with control `u`, we consider the **maximum** cost among all those `x⁺` states.
If the `cost_function` is `nothing`, we interpret this as a constant cost function of
`1` for each transition.

This function redirects to [`compute_worst_case_uniform_cost_controller`](@ref) if the
`cost_function` is nothing and [`compute_optimal_controller`](@ref) otherwise.
"""
function compute_worst_case_cost_controller(
    autom::ST.AbstractAutomatonList,
    target_set;
    initial_set = ST.enum_states(autom),
    cost_function = nothing,
    sparse_input::Bool = false,
)
    if cost_function === nothing
        return compute_worst_case_uniform_cost_controller(
            autom,
            target_set;
            initial_set = initial_set,
            sparse_input = sparse_input,
        )
    else
        return compute_optimal_controller(
            autom,
            target_set;
            initial_set = initial_set,
            sparse_input = sparse_input,
            cost_function = cost_function,
        )
    end
end

using DataStructures

"""
    compute_optimal_controller(
        autom::ST.AbstractAutomatonList,
        target_set;
        initial_set = ST.enum_states(autom),
        cost_function = nothing,
        sparse_input::Bool = false,
    )

Implementation for [`compute_worst_case_cost_controller`](@ref) supporting any
`cost_function` returning nonnegative values.
If the cost function returns the same value across all `x` and `u`, consider
calling the specialized function [`compute_worst_case_uniform_cost_controller`](@ref)
instead.
"""
function compute_optimal_controller(
    autom::ST.AbstractAutomatonList,
    target_set;
    initial_set = ST.enum_states(autom),
    cost_function = nothing,
    sparse_input::Bool = false,
)
    contr_tab = SymbolicControlTable(ST.get_n_state(autom))

    is_det = ST.is_deterministic(autom)
    uniform_cost = cost_function === nothing
    effective_cost_function = uniform_cost ? ((q, u) -> 1.0) : cost_function
    state_set = ST.enum_states(autom)
    value_fun_tab = fill(Inf, ST.get_n_state(autom))
    for q in target_set
        value_fun_tab[q] = 0.0
    end
    pq = PriorityQueue{Int, Float64}(q => 0.0 for q in target_set)

    num_init_unreachable = length(initial_set)
    optimal_controllable_set = Set(target_set)
    counter = is_det ? nothing : _counter(autom, sparse_input)

    while !isempty(pq) && num_init_unreachable > 0
        (target, cost_to_target) = dequeue_pair!(pq)

        if target in initial_set
            num_init_unreachable -= 1
        end
        if !is_det
            push!(optimal_controllable_set, target)
        end

        for (source, symbol) in ST.pre(autom, target)
            if is_det
                total_cost = effective_cost_function(source, symbol) + cost_to_target
                if total_cost < value_fun_tab[source]
                    value_fun_tab[source] = total_cost
                    add_control!(contr_tab, source, symbol)
                    pq[source] = total_cost
                end
            else
                if source in optimal_controllable_set
                    continue
                end
                if decrease_counter!(counter, source, symbol) == 0
                    if uniform_cost
                        worst = value_fun_tab[target]
                        push!(optimal_controllable_set, source)
                    else
                        successors = ST.post(autom, source, symbol)
                        worst = maximum(value_fun_tab[q′] for q′ in successors)
                    end
                    total_cost = effective_cost_function(source, symbol) + worst
                    if total_cost < value_fun_tab[source]
                        value_fun_tab[source] = total_cost
                        add_control!(contr_tab, source, symbol)
                        pq[source] = total_cost
                    end
                end
            end
        end
    end

    controllable_set = Set(i for (i, v) in pairs(value_fun_tab) if isfinite(v))
    uncontrollable_set = setdiff(state_set, controllable_set)
    controller = to_ms_controller(contr_tab)

    return controller, controllable_set, uncontrollable_set, value_fun_tab
end

function increase_counter!(counter::Array{Int, 2}, source::Int, symbol::Int)
    return counter[source, symbol] += 1
end
function increase_counter!(counter::Dict{Tuple{Int, Int}, Int}, source::Int, symbol::Int)
    key = (source, symbol)
    return counter[key] = get(counter, key, 0) + 1
end

function decrease_counter!(counter::Array{Int, 2}, source::Int, symbol::Int)
    counter[source, symbol] -= 1
    return counter[source, symbol]
end
function decrease_counter!(counter::Dict{Tuple{Int, Int}, Int}, source::Int, symbol::Int)
    key = (source, symbol)
    counter[key] = get(counter, key, 0) - 1
    return counter[key]
end

function _compute_num_targets_unreachable(counter, autom)
    for target in ST.enum_states(autom)
        for (source, symbol) in ST.pre(autom, target)
            increase_counter!(counter, source, symbol)
        end
    end
end

function _counter(autom, sparse_input::Bool)
    if sparse_input
        num_targets_unreachable = Dict{Tuple{Int, Int}, Int}()
    else
        num_targets_unreachable = zeros(Int, ST.get_n_state(autom), ST.get_n_input(autom))
    end

    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    return num_targets_unreachable
end

"""
    function compute_worst_case_uniform_cost_controller(
        autom::ST.AbstractAutomatonList,
        target_set;
        initial_set = ST.enum_states(autom),
        sparse_input = false,
    )

Implementation for [`compute_worst_case_uniform_cost_controller`](@ref) supporting for
a cost function that is `1` for any state `x` and input `u`.
Consider using [`compute_optimal_controller`](@ref) for a more general cost function.
But for this particular case of cost, this implementation is more efficient than [`compute_optimal_controller`](@ref)
as it does not rely on a `PriorityQueue`.
"""
function compute_worst_case_uniform_cost_controller(
    autom::ST.AbstractAutomatonList,
    target_set;
    initial_set = ST.enum_states(autom),
    sparse_input = false,
)
    contr_tab = SymbolicControlTable(ST.get_n_state(autom))

    stateset,
    initset,
    controllable_set,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab = _data(autom, initial_set, target_set, sparse_input)

    success, value_fun_tab = _compute_controller_reach!(
        contr_tab,
        autom,
        initset,
        controllable_set,
        num_targets_unreachable,
        current_targets,
        next_targets,
        value_fun_tab,
    )

    uncontrollable_set = setdiff(stateset, controllable_set)
    controller = to_ms_controller(contr_tab)

    return controller, controllable_set, uncontrollable_set, value_fun_tab
end

function _data(autom, initlist, targetlist, sparse_input::Bool)
    if sparse_input
        num_targets_unreachable = Dict{Tuple{Int, Int}, Int}()
    else
        num_targets_unreachable = zeros(Int, ST.get_n_state(autom), ST.get_n_input(autom))
    end

    _compute_num_targets_unreachable(num_targets_unreachable, autom)

    stateset = BitSet(ST.enum_states(autom))
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    value_fun_tab = fill(Inf, ST.get_n_state(autom)) # Inf = uncontrollable by default

    return stateset,
    initset,
    targetset,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab
end

function _compute_controller_reach!(
    contr_tab,
    autom,
    init_set,
    target_set,
    counter,
    current_targets,
    next_targets,
    value_fun_tab,
)::Tuple{Bool, Vector{Float64}}
    num_init_unreachable = length(init_set)

    step = 0
    for s in current_targets
        value_fun_tab[s] = step
    end

    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        step += 1

        for target in current_targets
            for (source, symbol) in ST.pre(autom, target)
                if !(source in target_set) &&
                   iszero(decrease_counter!(counter, source, symbol))
                    push!(target_set, source)
                    push!(next_targets, source)
                    add_control!(contr_tab, source, symbol)
                    value_fun_tab[source] = step

                    if source in init_set
                        num_init_unreachable -= 1
                    end
                end
            end
        end
        current_targets, next_targets = next_targets, current_targets
    end

    return iszero(num_init_unreachable), value_fun_tab
end
