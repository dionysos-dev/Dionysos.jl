
mutable struct OptimizerEmptyProblem{T} <: MOI.AbstractOptimizer
    # Input 
    empty_problem::Union{Nothing, PR.EmptyProblem}
    concrete_system::Union{Nothing, HybridSystem}
    optimizer_list::Union{Nothing, Any}
    optimizer_kwargs_dict::Union{Nothing, Any}

    max_iterations::Union{Nothing, Int}
    print_level::Int

    # Output 
    abstract_system::Union{Nothing, SY.TimedHybridSymbolicModel}
    abstraction_construction_time_sec::T

    function OptimizerEmptyProblem{T}() where {T}
        optimizer = new{T}(nothing, nothing, nothing, nothing, 1000, 1, nothing, 0.0)
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

function reset!(model::OptimizerEmptyProblem)
    return model
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

function MOI.optimize!(optimizer::OptimizerEmptyProblem)
    t_ref = time()
    # Ensure necessary parameters are set
    _validate_model(optimizer, [:empty_problem])
    @assert optimizer.empty_problem.system !== nothing "System must be set to construct the abstraction."

    # Build symbolic model
    optimizer.print_level >= 1 &&
        println("compute_abstract_system_from_concrete_system!: started")
    optimizer.abstract_system = SY.build_timed_hybrid_symbolic_model(
        optimizer.empty_problem.system,
        optimizer.optimizer_list,
        optimizer.optimizer_kwargs_dict,
    )

    optimizer.print_level >= 1 && println(
        "compute_abstract_system_from_concrete_system! terminated with success: ",
        "$(HybridSystems.ntransitions(SY.get_automaton(optimizer.abstract_system))) transitions created",
    )

    optimizer.abstraction_construction_time_sec = time() - t_ref
    return
end
