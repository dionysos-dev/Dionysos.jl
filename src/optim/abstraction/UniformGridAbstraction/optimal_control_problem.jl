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
    controllable_set::Union{Nothing, Dionysos.Domain.DomainList}
    uncontrollable_set::Union{Nothing, Dionysos.Domain.DomainList}

    function OptimizerOptimalControlProblem{T}() where {T}
        return new{T}(nothing, nothing, nothing, nothing, 0.0, true, nothing)
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
    return model.abstract_problem_time_sec  # ✅ Corrected: Ensures this tracks the actual solve time
end

function MOI.get(model::OptimizerOptimalControlProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

"""
    MOI.optimize!(optimizer::OptimizerOptimalControlProblem)

Solves the optimal control problem based on the abstract system.
"""
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
        Dionysos.Symbolic.enum_cells(optimizer.abstract_problem.system)

    abstract_controller, controllable_set_symbols, uncontrollable_set_symbols =
        compute_largest_controllable_set(
            optimizer.abstract_problem.system,
            optimizer.abstract_problem.target_set;
            initial_set = init_set,
        )

    controllable_set = Dionysos.Symbolic.get_domain_from_symbols(
        optimizer.abstract_system,
        controllable_set_symbols,
    )
    uncontrollable_set = Dionysos.Symbolic.get_domain_from_symbols(
        optimizer.abstract_system,
        uncontrollable_set_symbols,
    )

    optimizer.abstract_controller = abstract_controller
    optimizer.controllable_set = controllable_set
    optimizer.uncontrollable_set = uncontrollable_set

    # Display results
    if ⊆(optimizer.abstract_problem.initial_set, controllable_set_symbols)
        println("\n Reachability: terminated with success")
    else
        println("\n Reachability: terminated without covering initial set")
    end

    # Correctly track solve time
    return optimizer.abstract_problem_time_sec = time() - t_ref
end

"""
    build_abstract_problem(concrete_problem, abstract_system)

Builds the abstract problem based on the `concrete_problem` and `abstract_system`.
"""
function build_abstract_problem(
    concrete_problem::Dionysos.Problem.OptimalControlProblem,
    abstract_system::Dionysos.Symbolic.SymbolicModelList,
)
    @warn("The `state_cost` and `transition_cost` are not yet fully implemented")

    return Dionysos.Problem.OptimalControlProblem(
        abstract_system,
        Dionysos.Symbolic.get_symbols_from_set(
            abstract_system,
            concrete_problem.initial_set,
            Dionysos.Domain.OUTER,
        ),
        Dionysos.Symbolic.get_symbols_from_set(
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
)
    abstract_controller = NewControllerList()
    stateset,
    initset,
    controllable_set,
    num_targets_unreachable,
    current_targets,
    next_targets = _data(abstract_system.autom, initial_set, target_set)

    _compute_controller_reach!(
        abstract_controller,
        abstract_system.autom,
        initset,
        controllable_set,
        num_targets_unreachable,
        current_targets,
        next_targets,
    )

    uncontrollable_set = setdiff(stateset, controllable_set)

    return abstract_controller, controllable_set, uncontrollable_set
end

"""
    _compute_num_targets_unreachable(num_targets_unreachable, autom)

Computes the number of targets unreachable based on the `autom`.
"""
function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    for target in 1:(autom.nstates)
        for soursymb in Dionysos.Symbolic.pre(autom, target)
            num_targets_unreachable[soursymb[1], soursymb[2]] += 1
        end
    end
end

"""
    _data(contr, autom, initlist, targetlist)

Returns the data based on the `autom`, `initlist`, and `targetlist`.
"""
function _data(autom, initlist, targetlist)
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    stateset = BitSet(1:(autom.nstates))
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    return stateset,
    initset,
    targetset,
    num_targets_unreachable,
    current_targets,
    next_targets
end

"""
    _compute_controller_reach!(contr, autom, init_set, target_set, num_targets_unreachable, current_targets, next_targets)

Computes the controller reach based on the `contr`, `autom`, `init_set`, `target_set`, `num_targets_unreachable`, `current_targets`, and `next_targets`.
"""
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
