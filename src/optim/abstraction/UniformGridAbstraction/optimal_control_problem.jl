mutable struct OptimizerOptimalControlProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem}
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}

    # Common parameters
    abstract_problem::Union{Nothing, Dionysos.Problem.OptimalControlProblem}
    abstract_controller::Union{Nothing, Dionysos.Utils.SortedTupleSet{2, NTuple{2, Int}}}
    abstract_problem_time_sec::T

    # Specific parameters
    early_stop::Union{Nothing, Bool}
    sparse_input::Bool
    controllable_set::Union{Nothing, Dionysos.Domain.DomainList}
    uncontrollable_set::Union{Nothing, Dionysos.Domain.DomainList}
    value_fun_tab::Union{Nothing, Any} # Value function in tabular form, Inf means uncontrollable state
    abstract_value_function::Union{Nothing, Any}
    concrete_value_function::Union{Nothing, Any}

    success::Bool
    print_level::Int

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            0.0,
            true,
            false,
            nothing,
            nothing,
            false,
            1,
        )
    end
end

OptimizerOptimalControlProblem() = OptimizerOptimalControlProblem{Float64}()

MOI.is_empty(optimizer::OptimizerOptimalControlProblem) =
    optimizer.concrete_problem === nothing

function MOI.set(
    model::OptimizerOptimalControlProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerOptimalControlProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerOptimalControlProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstract_value_function(value_fun_tab)
    return abstract_value_function(state) = value_fun_tab[state]
end

function build_concrete_value_function(abstract_system, abstract_value_function)
    function concrete_value_function(x)
        state = SY.get_abstract_state(abstract_system, x)
        return abstract_value_function(state)
    end
    return concrete_value_function
end

function MOI.optimize!(optimizer::OptimizerOptimalControlProblem)
    t_ref = time()

    # Ensure abstract system is set before proceeding
    if optimizer.abstract_system === nothing
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    end

    # Build the abstract problem from the concrete problem
    optimizer.abstract_problem =
        build_abstract_problem(optimizer.concrete_problem, optimizer.abstract_system)

    # Compute the largest controllable set
    init_set =
        optimizer.early_stop ? optimizer.abstract_problem.initial_set :
        Dionysos.Symbolic.enum_states(optimizer.abstract_problem.system)

    optimizer.print_level >= 1 && println("compute_controller_reachability! started")
    abstract_controller,
    controllable_set_symbols,
    uncontrollable_set_symbols,
    value_fun_tab = compute_largest_controllable_set(
        optimizer.abstract_problem.system,
        optimizer.abstract_problem.target_set;
        initial_set = init_set,
        sparse_input = optimizer.sparse_input,
    )

    controllable_set = Dionysos.Symbolic.get_domain_from_states(
        optimizer.abstract_system,
        controllable_set_symbols,
    )
    uncontrollable_set = Dionysos.Symbolic.get_domain_from_states(
        optimizer.abstract_system,
        uncontrollable_set_symbols,
    )

    optimizer.abstract_controller = abstract_controller
    optimizer.controllable_set = controllable_set
    optimizer.uncontrollable_set = uncontrollable_set
    optimizer.value_fun_tab = value_fun_tab
    optimizer.abstract_value_function =
        build_abstract_value_function(optimizer.value_fun_tab)
    optimizer.concrete_value_function = build_concrete_value_function(
        optimizer.abstract_system,
        optimizer.abstract_value_function,
    )

    # Display results
    if ⊆(optimizer.abstract_problem.initial_set, controllable_set_symbols)
        optimizer.success = true
    end
    optimizer.print_level >= 1 &&
        println("\n Reachability: terminated with $(optimizer.success)")
    optimizer.abstract_problem_time_sec = time() - t_ref
    return
end

function build_abstract_problem(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
)
    @warn("The `state_cost` and `transition_cost` are not yet fully implemented")

    return Dionysos.Problem.OptimalControlProblem(
        abstract_system,
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.initial_set,
            Dionysos.Domain.OUTER,
        ),
        Dionysos.Symbolic.get_states_from_set(
            abstract_system,
            concrete_problem.target_set,
            Dionysos.Domain.INNER,
        ),
        concrete_problem.state_cost,       # TODO: Transform continuous cost into discrete abstraction
        concrete_problem.transition_cost,  # TODO: Transform continuous cost into discrete abstraction
        concrete_problem.time,              # TODO: Translate continuous time into discrete steps
    )
end

"""
    compute_largest_controllable_set(abstract_system, target_set; initial_set)

Computes the largest controllable set and generates an abstract controller.
"""
function compute_largest_controllable_set(
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
    target_set;
    initial_set = Dionysos.Symbolic.enum_cells(abstract_system),
    sparse_input = false,
)
    abstract_controller = NewControllerList()
    stateset,
    initset,
    controllable_set,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab = _data(abstract_system.autom, initial_set, target_set, sparse_input)

    success, value_fun_tab = _compute_controller_reach!(
        abstract_controller,
        abstract_system.autom,
        initset,
        controllable_set,
        num_targets_unreachable,
        current_targets,
        next_targets,
        value_fun_tab,
    )

    uncontrollable_set = setdiff(stateset, controllable_set)

    return abstract_controller, controllable_set, uncontrollable_set, value_fun_tab
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
    for target in 1:(autom.nstates)
        for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
            increase_counter!(counter, source, symbol)
        end
    end
end

function _data(autom, initlist, targetlist, sparse_input::Bool)
    if sparse_input
        num_targets_unreachable = Dict{Tuple{Int, Int}, Int}()
    else
        num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols)
    end

    _compute_num_targets_unreachable(num_targets_unreachable, autom)

    stateset = BitSet(1:(autom.nstates))
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    value_fun_tab = fill(Inf, autom.nstates) # Inf = uncontrollable by default

    return stateset,
    initset,
    targetset,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab
end

function _compute_controller_reach!(
    contr,
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
            for (source, symbol) in Dionysos.Symbolic.pre(autom, target)
                if !(source in target_set) &&
                   iszero(decrease_counter!(counter, source, symbol))
                    push!(target_set, source)
                    push!(next_targets, source)
                    Dionysos.Utils.push_new!(contr, (source, symbol))
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
