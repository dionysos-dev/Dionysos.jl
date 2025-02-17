@enum ApproxMode GROWTH LINEARIZED DELTA_GAS

mutable struct OptimizerEmptyProblem{T} <: MOI.AbstractOptimizer
    # Metaparameters
    state_grid::Union{Nothing, Dionysos.Domain.Grid}
    input_grid::Union{Nothing, Dionysos.Domain.Grid}
    empty_problem::Union{Nothing, Dionysos.Problem.EmptyProblem}

    # Constructed meta-parameters
    discretized_system::Any
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}
    abstraction_construction_time_sec::T

    # Reachable_function
    approx_mode::ApproxMode
    δGAS::Union{Nothing, Bool}

    ## Discrete time system
    growthbound_map::Union{Nothing, Function}
    sys_inv_map::Union{Nothing, Function}

    ## Continuous time system
    time_step::T
    num_sub_steps_system_map::Int
    num_sub_steps_growth_bound::Int

    # Specific reachable_function based on growth bound
    jacobian_bound::Union{Nothing, Function}

    function OptimizerEmptyProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0.0,
            GROWTH,
            false,
            nothing,
            nothing,
            NaN,
            5,
            5,
            nothing,
        )
    end
end

OptimizerEmptyProblem() = OptimizerEmptyProblem{Float64}()

MOI.is_empty(optimizer::OptimizerEmptyProblem) = optimizer.empty_problem === nothing

function MOI.set(model::OptimizerEmptyProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerEmptyProblem, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec  # ✅ Returns abstraction construction time
end

function MOI.get(model::OptimizerEmptyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

"""
    validate_model!(model::OptimizerEmptyProblem, required_fields::Vector{Symbol})

Ensures the specified `required_fields` are set in the `model`. Throws an error
if any field is missing.
"""
function _validate_model(model::OptimizerEmptyProblem, required_fields::Vector{Symbol})
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set the `$(field)`. Missing required field in OptimizerEmptyProblem.",
            )
        end
    end
end

"""
    validate_discrete_model!(model::OptimizerEmptyProblem)

Validates that the `OptimizerEmptyProblem` model is correctly configured for discrete
systems. Throws an error if required fields like `sys_inv_map` or `growthbound_map`
are missing or invalid.
"""
function _validate_discrete_model(model::OptimizerEmptyProblem)
    return _validate_model(model, [:growthbound_map, :sys_inv_map])
end

"""
    validate_continuous_model!(model::OptimizerEmptyProblem)

Validates that the `OptimizerEmptyProblem` model is correctly configured for continuous
systems. Throws an error if required fields like `time_step` or `jacobian_bound`
are missing or invalid.
"""
function _validate_continuous_model(model::OptimizerEmptyProblem)
    _validate_model(model, [:time_step, :jacobian_bound])
    if isnan(model.time_step)
        error("Please set a valid `time_step`.")
    end
end

"""
    _maybe_discretized_system(concrete_system, model, noise)

Returns the discretized system based on the `concrete_system` and `model`.
"""
function _maybe_discretized_system(concrete_system, model, noise)
    if isa(concrete_system, MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem)
        _validate_discrete_model(model)
        return _discrete_system(concrete_system, model, noise)
    elseif isa(
        concrete_system,
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem,
    )
        _validate_continuous_model(model)
        return _discretize_continuous_system(concrete_system, model, noise)
    else
        error("Unsupported system type: $(typeof(concrete_system))")
    end
end

"""
    discretized_system(concrete_system::MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem, model, noise)

Returns the discretized system based on the `concrete_system` and `model`.
"""
function _discrete_system(
    concrete_system::MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem,
    model,
    noise,
)
    return Dionysos.System.ControlSystemGrowth(
        1.0, # `time_step` should be ignored by `concrete_system.f`, `model.growthbound_map` and `model.sys_inv_map` anyway
        noise,
        noise,
        concrete_system.f,
        model.growthbound_map,
        model.sys_inv_map,
    )
end

"""
    _discretize_continuous_system(concrete_system::MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem, model, noise)

Returns the discretized system based on the `concrete_system` and `model`.
"""
function _discretize_continuous_system(
    concrete_system::MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem,
    model,
    noise,
)
    if model.approx_mode == GROWTH
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

"""
    _get_domain_list(model, variables, grid)

Returns a `Dionysos.Domain.DomainList` based on the `variables` and `grid`.
"""
function _get_domain_list(variables, grid, position_in_domain)
    domain_list = Dionysos.Domain.DomainList(grid)
    Dionysos.Domain.add_set!(domain_list, variables, position_in_domain)
    return domain_list
end

_vector_of_tuple(size, value = 0.0) = StaticArrays.SVector(ntuple(_ -> value, Val(size)))

function MOI.optimize!(optimizer::OptimizerEmptyProblem)
    t_ref = time()

    # Ensure necessary parameters are set
    _validate_model(optimizer, [:state_grid, :input_grid, :empty_problem])

    concrete_system = optimizer.empty_problem.system

    # Create abstract system
    abstract_system = Dionysos.Symbolic.NewSymbolicModelListList(
        _get_domain_list(concrete_system.X, optimizer.state_grid, Dionysos.Domain.INNER),
        _get_domain_list(concrete_system.U, optimizer.input_grid, Dionysos.Domain.CENTER),
    )

    # TODO: Consider adding noise handling
    @warn("Noise is not yet accounted for in system abstraction.")
    noise = _vector_of_tuple(Dionysos.Utils.get_dims(concrete_system.X))

    optimizer.discretized_system =
        _maybe_discretized_system(concrete_system, optimizer, noise)

    if optimizer.δGAS
        Dionysos.Symbolic.compute_deterministic_symmodel_from_controlsystem!(
            abstract_system,
            optimizer.discretized_system,
        )
    else
        Dionysos.Symbolic.compute_symmodel_from_controlsystem!(
            abstract_system,
            optimizer.discretized_system,
        )
    end

    optimizer.abstract_system = abstract_system
    return optimizer.abstraction_construction_time_sec = time() - t_ref  # Track the abstraction time
end

##### Mulitdipatch for compute_reachable_set.

"""
    get_compute_reachable_set(concrete_system) -> Function

Determines the appropriate `compute_reachable_set` function based on the provided system data.
- If the user explicitly provides one, use it.
- If `concrete_system` includes a `growthbound_map`, use it.
- If `concrete_system` includes `jacobian_bound`, compute the `growthbound_map` first.
- If only dynamics are provided, estimate Lipschitz bounds and derive `growthbound_map`.
"""
function get_compute_reachable_set(concrete_system)
    # 1️⃣ If the system already provides `compute_reachable_set`, use it.
    if hasmethod(concrete_system, (:compute_reachable_set,))
        return concrete_system.compute_reachable_set
    end

    # 2️⃣ If the system is discrete-time and has `growthbound_map`, use it directly.
    if isa(concrete_system, ST.ControlSystemGrowth)
        return (rect, u) -> compute_reachable_set(concrete_system, rect, u)
    end

    # 3️⃣ If the system has a `jacobian_bound`, derive `growthbound_map` first.
    if hasproperty(concrete_system, :jacobian_bound)
        derived_growthbound_map = derive_growth_bound(concrete_system)
        return (rect, u) -> compute_reachable_set_with_growth_bound(
            concrete_system,
            rect,
            u,
            derived_growthbound_map,
        )
    end

    # 4️⃣ If the system only provides dynamics, estimate Lipschitz bounds first.
    if hasproperty(concrete_system, :dynamics)
        estimated_growthbound_map = estimate_growth_bound(concrete_system)
        return (rect, u) -> compute_reachable_set_with_growth_bound(
            concrete_system,
            rect,
            u,
            estimated_growthbound_map,
        )
    end

    return error(
        "Cannot determine a valid compute_reachable_set function for the given system.",
    )
end

"""
    compute_reachable_set(concrete_system, rect, u)

Computes an over-approximation of the reachable set.
"""
function compute_reachable_set(
    concrete_system::ST.ControlSystemGrowth{N},
    rect::UT.HyperRectangle,
    u,
) where {N}
    x = UT.get_center(rect)
    r = UT.get_r(rect)
    Fx = concrete_system.sys_map(x, u, concrete_system.tstep)
    Fr =
        concrete_system.growthbound_map(r, u, concrete_system.tstep) +
        concrete_system.measnoise
    return UT.HyperRectangle(Fx - Fr, Fx + Fr)
end

"""
    derive_growth_bound(concrete_system) -> Function

Derives a growth bound function from the Jacobian.
"""
function derive_growth_bound(concrete_system)
    return (r, u, t) -> begin
        J = concrete_system.jacobian_bound(u)
        return exp(t * J) * r
    end
end

"""
    estimate_growth_bound(concrete_system) -> Function

Estimates a growth bound function from the Lipschitz constant.
"""
function estimate_growth_bound(concrete_system)
    L = estimate_lipschitz_constant(concrete_system)
    return (r, u, t) -> begin
        return exp(t * L) * r
    end
end

"""
    estimate_lipschitz_constant(concrete_system) -> Real

Computes an estimate of the Lipschitz constant.
"""
function estimate_lipschitz_constant(concrete_system)
    # Approximate Lipschitz constant using finite differences or heuristics.
    return maximum(abs, concrete_system.dynamics)
end

# Feature	Solution
# User provides a custom compute_reachable_set	We use it directly.
# System has a growthbound_map	Use it for reachable set computation.
# System has jacobian_bound	Compute growthbound_map first, then use it.
# Only system dynamics are available	Estimate Lipschitz bound, then use it.
# System is continuous-time	Convert it to discrete-time first.

# ✔ Modular → The system dynamically chooses the best method.
# ✔ Efficient → Uses the most accurate information available.
