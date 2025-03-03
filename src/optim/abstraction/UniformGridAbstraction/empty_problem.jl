@enum ApproxMode USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION

mutable struct OptimizerEmptyProblem{T} <: MOI.AbstractOptimizer
    ### Abstraction Result
    discrete_time_system::Union{
        Nothing,
        MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem,
    }
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}
    abstraction_construction_time_sec::T

    ### System Approximation
    continuous_time_system_approximation::Union{
        Nothing,
        ST.ContinuousTimeSystemApproximation,
    }
    discrete_time_system_approximation::Union{Nothing, ST.DiscreteTimeSystemApproximation}

    ### User Settings
    empty_problem::Union{Nothing, Dionysos.Problem.EmptyProblem}
    state_grid::Union{Nothing, Dionysos.Domain.Grid}
    input_grid::Union{Nothing, Dionysos.Domain.Grid}

    ### USER_DEFINED
    overapproximation_map::Union{Nothing, Function}

    ### GROWTH
    growthbound_map::Union{Nothing, Function}
    jacobian_bound::Union{Nothing, Function}
    ngrowthbound::Int

    ### LINEARIZED
    DF_sys::Union{Nothing, Function}  # Jacobian function
    bound_DF::Union{Nothing, Function}  # Bound on Jacobian
    bound_DDF::Union{Nothing, Function}  # Bound on Hessian

    ### RANDOM_SIMULATION
    n_samples::Union{Nothing, Int}

    ### Continuous-Time System Settings
    time_step::Union{Nothing, T}
    nsystem::Int

    ### System Approximation Settings
    approx_mode::ApproxMode
    efficient::Bool

    print_level::Int

    function OptimizerEmptyProblem{T}() where {T}
        optimizer = new{T}(
            nothing,
            nothing,
            0.0, # Abstraction
            nothing,
            nothing, # System Approximation
            nothing,
            nothing,
            nothing, # User Settings
            nothing, # USER_DEFINED
            nothing,
            nothing,
            5, # GROWTH
            nothing,
            nothing,
            nothing, # LINEARIZED
            nothing,  # RANDOM_SIMULATION
            nothing,
            5,  # Continuous-Time System
            GROWTH,
            true, # Approximation Settings
            1,
        )
        return optimizer
    end
end

OptimizerEmptyProblem() = OptimizerEmptyProblem{Float64}()

MOI.is_empty(optimizer::OptimizerEmptyProblem) = optimizer.empty_problem === nothing

function MOI.set(model::OptimizerEmptyProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerEmptyProblem, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec
end

function MOI.get(model::OptimizerEmptyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function _validate_model(model::OptimizerEmptyProblem, required_fields::Vector{Symbol})
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set the `$(field)`. Missing required field in OptimizerEmptyProblem.",
            )
        end
    end
end

function build_continuous_approximation(optimizer::OptimizerEmptyProblem, system)
    mode = optimizer.approx_mode
    if mode == USER_DEFINED
        _validate_model(optimizer, [:overapproximation_map])
        return ST.ContinuousTimeOverapproximationMap(
            system,
            optimizer.overapproximation_map,
        )
    elseif mode == GROWTH
        if optimizer.growthbound_map !== nothing
            return ST.ContinuousTimeGrowthBound(system, optimizer.growthbound_map)
        elseif optimizer.jacobian_bound !== nothing
            return ST.ContinuousTimeGrowthBound_from_jacobian_bound(
                system,
                optimizer.jacobian_bound,
            )
        else
            return ST.ContinuousTimeGrowthBound(system)
        end
    elseif mode == LINEARIZED
        _validate_model(optimizer, [:DF_sys, :bound_DF, :bound_DDF])
        return ST.ContinuousTimeLinearized(
            system,
            optimizer.DF_sys,
            optimizer.bound_DF,
            optimizer.bound_DDF,
        )
    elseif mode == CENTER_SIMULATION
        return ST.ContinuousTimeCenteredSimulation(system)
    elseif mode == RANDOM_SIMULATION
        _validate_model(optimizer, [:n_samples])
        return ST.ContinuousTimeRandomSimulation(system, optimizer.n_samples)

    else
        error("Unsupported approximation mode: $mode")
    end
end

function build_discrete_approximation(optimizer::OptimizerEmptyProblem, system)
    mode = optimizer.approx_mode
    if mode == USER_DEFINED
        _validate_model(optimizer, [:overapproximation_map])
        return ST.DiscreteTimeOverapproximationMap(system, optimizer.overapproximation_map)
    elseif mode == GROWTH
        _validate_model(optimizer, [:growthbound_map])
        return ST.DiscreteTimeGrowthBound(system, optimizer.growthbound_map)
    elseif mode == LINEARIZED
        _validate_model(optimizer, [:DF_sys, :bound_DF, :bound_DDF])
        return ST.DiscreteTimeLinearized(
            system,
            optimizer.DF_sys,
            optimizer.bound_DF,
            optimizer.bound_DDF,
        )
    elseif mode == CENTER_SIMULATION
        return ST.DiscreteTimeCenteredSimulation(system)
    elseif mode == RANDOM_SIMULATION
        _validate_model(optimizer, [:n_samples])
        return ST.DiscreteTimeRandomSimulation(system, optimizer.n_samples)
    else
        error("Unsupported approximation mode: $mode")
    end
end

function build_system_approximation!(optimizer::OptimizerEmptyProblem)
    _validate_model(optimizer, [:empty_problem])
    @assert optimizer.empty_problem.system !== nothing "System must be set before building overapproximation."

    system = optimizer.empty_problem.system
    if isa(system, MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem)
        _validate_model(optimizer, [:time_step])  # Ensure time step is provided
        optimizer.continuous_time_system_approximation =
            build_continuous_approximation(optimizer, system)
        optimizer.discrete_time_system_approximation = ST.discretize(
            optimizer.continuous_time_system_approximation,
            optimizer.time_step,
        )
    elseif isa(system, MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem)
        optimizer.discrete_time_system_approximation =
            build_discrete_approximation(optimizer, system)
    else
        error("Unknown system type: $(typeof(system))")
    end
    return optimizer.discrete_time_system =
        ST.get_system(optimizer.discrete_time_system_approximation)
end

function _get_domain_list(variables, grid, position_in_domain)
    domain_list = Dionysos.Domain.DomainList(grid)
    Dionysos.Domain.add_set!(domain_list, variables, position_in_domain)
    return domain_list
end

_vector_of_tuple(size, value = 0.0) = SVector(ntuple(_ -> value, Val(size)))

function MOI.optimize!(optimizer::OptimizerEmptyProblem)
    t_ref = time()

    # Ensure necessary parameters are set
    _validate_model(optimizer, [:state_grid, :input_grid, :empty_problem])
    @assert optimizer.empty_problem.system !== nothing "System must be set before building overapproximation."

    # Create over-approximation method
    build_system_approximation!(optimizer)

    concrete_system = optimizer.empty_problem.system

    # Create abstract system
    abstract_system = Dionysos.Symbolic.NewSymbolicModelListList(
        _get_domain_list(concrete_system.X, optimizer.state_grid, Dionysos.Domain.INNER),
        _get_domain_list(concrete_system.U, optimizer.input_grid, Dionysos.Domain.CENTER),
    )

    if optimizer.print_level >= 2
        @info("Number of states: $(SY.get_n_state(abstract_system))")
        @info("Number of inputs: $(SY.get_n_input(abstract_system))")
        @info(
            "Number of forward images: $(SY.get_n_input(abstract_system)*SY.get_n_state(abstract_system))"
        )
    end

    # TODO: Consider adding noise handling
    @warn("Noise is not yet accounted for in system abstraction.")
    noise = _vector_of_tuple(Dionysos.Utils.get_dims(concrete_system.X))

    optimizer.print_level >= 1 && println(
        "compute_abstract_system_from_concrete_system!: started with $(typeof(optimizer.discrete_time_system_approximation))",
    )
    if !optimizer.efficient &&
       ST.is_over_approximation(optimizer.discrete_time_system_approximation)
        Dionysos.Symbolic.compute_abstract_system_from_concrete_system!(
            abstract_system,
            ST.get_DiscreteTimeOverApproximationMap(
                optimizer.discrete_time_system_approximation,
            );
            verbose = optimizer.print_level >= 2,
        )
    else
        Dionysos.Symbolic.compute_abstract_system_from_concrete_system!(
            abstract_system,
            optimizer.discrete_time_system_approximation;
            verbose = optimizer.print_level >= 2,
        )
    end
    optimizer.print_level >= 1 && println(
        "compute_abstract_system_from_concrete_system! terminated with success: ",
        "$(HybridSystems.ntransitions(abstract_system.autom)) transitions created",
    )

    optimizer.abstract_system = abstract_system
    optimizer.abstraction_construction_time_sec = time() - t_ref
    return
end
