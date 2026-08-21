# ============================================================
# Optimal Control
# ============================================================

"""
    OptimizerOptimalControlProblem{T} <: AbstractDionysosOptimizer

Reach-avoid synthesis on a finite automaton: given an
[`OptimalControlProblem`](@ref Dionysos.Problem.OptimalControlProblem) whose system is an
`AbstractAutomatonList`, compute the controllable set and a controller reaching the target
without leaving the safe set.

Which algorithm runs depends on the cost: unit cost sweeps backward in breadth-first layers, a
general nonnegative cost uses a priority queue. Setting `"bounded_input_variation"` instead
routes the synthesis to [`compute_bounded_input_variation_controller`](@ref), which constrains
*consecutive* inputs and returns a dynamic controller.

Set `"problem"`, optionally `"early_stop"`, `"sparse_input"`, `"bounded_input_variation"`; read
back `"controller"`, `"controllable_set"`, `"uncontrollable_set"` and `"value_function"`.
"""
mutable struct OptimizerOptimalControlProblem{T} <: AbstractDionysosOptimizer
    # inputs
    problem::Union{Nothing, PR.OptimalControlProblem}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int
    # Optional input slew-rate constraint; routes the synthesis to
    # `compute_bounded_input_variation_controller`.
    bounded_input_variation::Union{Nothing, BoundedInputVariation}

    # outputs
    controller::Union{Nothing, ST.AbstractDiscreteController}
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
            nothing, # bounded_input_variation
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

function build_value_function(value_fun_tab)
    return state -> value_fun_tab[state]
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system
    init_set = optimizer.early_stop ? problem.initial_set : SY.enum_states(autom)

    if optimizer.bounded_input_variation !== nothing
        controller, controllable_set, uncontrollable_set, value_fun_tab =
            compute_bounded_input_variation_controller(
                autom,
                problem.target_set,
                optimizer.bounded_input_variation;
                initial_set = init_set,
                safe_set = problem.safe_set,
                cost_function = problem.transition_cost,
            )
    else
        controller, controllable_set, uncontrollable_set, value_fun_tab =
            compute_worst_case_cost_controller(
                autom,
                problem.target_set;
                initial_set = init_set,
                safe_set = problem.safe_set,
                sparse_input = optimizer.sparse_input,
                cost_function = problem.transition_cost,
            )
    end

    optimizer.controller = controller
    optimizer.controllable_set = controllable_set
    optimizer.uncontrollable_set = uncontrollable_set
    optimizer.value_fun_tab = value_fun_tab
    optimizer.value_function = build_value_function(value_fun_tab)
    optimizer.success = covers_initial_set(q -> q in controllable_set, problem.initial_set)

    optimizer.print_level >= 1 &&
        println("Optimal control terminated with success = $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end

"""
    compute_worst_case_cost_controller(
        autom::SY.AbstractAutomatonList,
        target_set;
        initial_set = SY.enum_states(autom),
        safe_set = nothing,
        cost_function = nothing,
        sparse_input::Bool = false,
    )

Compute controller for optimal control problem assuming **worst case** for uncertainty.
This means, if there are several states `x⁺` that can be the next state from state
`x` with control `u`, we consider the **maximum** cost among all those `x⁺` states.
If the `cost_function` is `nothing`, we interpret this as a constant cost function of
`1` for each transition.

`safe_set` restricts the synthesis to a reach-avoid specification: only states in it may
enter the controllable set. `nothing` means the whole state space.

This function redirects to [`compute_worst_case_uniform_cost_controller`](@ref) if the
`cost_function` is nothing and [`compute_optimal_controller`](@ref) otherwise.
"""
function compute_worst_case_cost_controller(
    autom::SY.AbstractAutomatonList,
    target_set;
    initial_set = SY.enum_states(autom),
    safe_set = nothing,
    cost_function = nothing,
    sparse_input::Bool = false,
)
    if cost_function === nothing
        return compute_worst_case_uniform_cost_controller(
            autom,
            target_set;
            initial_set = initial_set,
            safe_set = safe_set,
            sparse_input = sparse_input,
        )
    else
        return compute_optimal_controller(
            autom,
            target_set;
            initial_set = initial_set,
            safe_set = safe_set,
            sparse_input = sparse_input,
            cost_function = cost_function,
        )
    end
end

# The safe set as a dense bit mask, so the synthesis loops stay type-stable and branch on a
# plain `Bool`. `nothing` — no avoid part — becomes an all-true mask rather than a second code
# path, which is why the fixed points below need no `safe_set === nothing` special case.
_safe_bits(safe_set, nstates::Int) =
    safe_set === nothing ? trues(nstates) : _bitset_from_states(safe_set, nstates)

# A target state outside the safe set does not count as reached: the specification is
# `safe U target`, so the safe part must still hold when the target is entered.
_safe_targets(target_set, safe_bits::BitVector) =
    Iterators.filter(q -> safe_bits[q], target_set)

using DataStructures

"""
    compute_optimal_controller(
        autom::SY.AbstractAutomatonList,
        target_set;
        initial_set = SY.enum_states(autom),
        safe_set = nothing,
        cost_function = nothing,
        sparse_input::Bool = false,
    )

Implementation for [`compute_worst_case_cost_controller`](@ref) supporting any
`cost_function` returning nonnegative values. `safe_set` is the optional avoid part of a
reach-avoid specification; see [`compute_worst_case_cost_controller`](@ref).
If the cost function returns the same value across all `x` and `u`, consider
calling the specialized function [`compute_worst_case_uniform_cost_controller`](@ref)
instead.
"""
function compute_optimal_controller(
    autom::SY.AbstractAutomatonList,
    target_set;
    initial_set = SY.enum_states(autom),
    safe_set = nothing,
    cost_function = nothing,
    sparse_input::Bool = false,
)
    contr_tab = DiscreteControlTable(SY.get_n_state(autom))

    is_det = SY.is_deterministic(autom)
    uniform_cost = cost_function === nothing
    effective_cost_function = uniform_cost ? ((q, u) -> 1.0) : cost_function
    state_set = SY.enum_states(autom)
    safe_bits = _safe_bits(safe_set, SY.get_n_state(autom))
    seeds = collect(_safe_targets(target_set, safe_bits))
    value_fun_tab = fill(Inf, SY.get_n_state(autom))
    for q in seeds
        value_fun_tab[q] = 0.0
    end
    pq = PriorityQueue{Int, Float64}(q => 0.0 for q in seeds)

    num_init_unreachable = length(initial_set)
    optimal_controllable_set = Set(seeds)
    counter = is_det ? nothing : _counter(autom, sparse_input)

    while !isempty(pq) && num_init_unreachable > 0
        (target, cost_to_target) = dequeue_pair!(pq)

        if target in initial_set
            num_init_unreachable -= 1
        end
        if !is_det
            push!(optimal_controllable_set, target)
        end

        for (source, symbol) in SY.pre(autom, target)
            # An unsafe state can never be part of a run that satisfies `safe U target`, so it
            # never enters the controllable set. Any input that may land in it therefore never
            # sees its counter reach zero, and is never selected — which is precisely the
            # worst-case avoidance the specification asks for.
            safe_bits[source] || continue
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
                        successors = SY.post(autom, source, symbol)
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
    controller = ST.DiscreteStaticController(controllable_set, contr_tab, false)

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
    for target in SY.enum_states(autom)
        for (source, symbol) in SY.pre(autom, target)
            increase_counter!(counter, source, symbol)
        end
    end
end

function _counter(autom, sparse_input::Bool)
    if sparse_input
        num_targets_unreachable = Dict{Tuple{Int, Int}, Int}()
    else
        num_targets_unreachable = zeros(Int, SY.get_n_state(autom), SY.get_n_input(autom))
    end

    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    return num_targets_unreachable
end

"""
    function compute_worst_case_uniform_cost_controller(
        autom::SY.AbstractAutomatonList,
        target_set;
        initial_set = SY.enum_states(autom),
        safe_set = nothing,
        sparse_input = false,
    )

Implementation for [`compute_worst_case_uniform_cost_controller`](@ref) supporting for
a cost function that is `1` for any state `x` and input `u`.
Consider using [`compute_optimal_controller`](@ref) for a more general cost function.
But for this particular case of cost, this implementation is more efficient than [`compute_optimal_controller`](@ref)
as it does not rely on a `PriorityQueue`.
"""
function compute_worst_case_uniform_cost_controller(
    autom::SY.AbstractAutomatonList,
    target_set;
    initial_set = SY.enum_states(autom),
    safe_set = nothing,
    sparse_input = false,
)
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    contr_tab = DiscreteControlTable(nstates)

    safe_bits = _safe_bits(safe_set, nstates)
    init_bits = _bitset_from_states(initial_set, nstates)
    controllable_bits = _bitset_from_states(target_set, nstates)
    controllable_bits .&= safe_bits

    counter = sparse_input ? _counter(autom, true) : _counter_dense(autom)

    current_targets = Int[]
    for q in _safe_targets(target_set, safe_bits)
        push!(current_targets, q)
    end

    next_targets = Int[]
    value_fun_tab = fill(Inf, nstates)

    step = 0
    for q in current_targets
        value_fun_tab[q] = 0.0
    end

    num_init_unreachable = 0
    @inbounds for q in 1:nstates
        init_bits[q] && !controllable_bits[q] && (num_init_unreachable += 1)
    end

    while !isempty(current_targets) && num_init_unreachable > 0
        empty!(next_targets)
        step += 1

        for target in current_targets
            for (source, symbol) in SY.pre(autom, target)
                controllable_bits[source] && continue
                # See `compute_optimal_controller`: unsafe states stay out of the controllable
                # set, which blocks every input that risks reaching them.
                safe_bits[source] || continue

                if decrease_counter!(counter, source, symbol) == 0
                    controllable_bits[source] = true
                    push!(next_targets, source)

                    add_control!(contr_tab, source, symbol)
                    value_fun_tab[source] = step

                    if init_bits[source]
                        num_init_unreachable -= 1
                    end
                end
            end
        end

        current_targets, next_targets = next_targets, current_targets
    end

    controllable_set = _set_from_bitset(controllable_bits)

    uncontrollable_bits = .!controllable_bits
    uncontrollable_set = _set_from_bitset(uncontrollable_bits)

    controller = ST.DiscreteStaticController(controllable_set, contr_tab, false)

    return controller, controllable_set, uncontrollable_set, value_fun_tab
end
