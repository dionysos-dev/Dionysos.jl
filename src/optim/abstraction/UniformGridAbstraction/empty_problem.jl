@enum ApproxMode USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION

"""
    OptimizerEmptyProblem{T} <: MOI.AbstractOptimizer

A solver responsible for constructing an **abstraction of the system dynamics**, independently of the control specification.

This optimizer wraps everything needed to solve an [`EmptyProblem`](@ref Dionysos.Problem.EmptyProblem), which is used to generate a symbolic model (abstraction) of either a continuous- or discrete-time system.

---

### Purpose

This optimizer builds a symbolic abstraction by simulating or approximating the behavior of the given system.  
The abstraction method is chosen via the `approx_mode` field, and determines which parameters and approximation logic are used.

---

### Parameters

#### Mandatory fields set by the user

- `empty_problem` (**required**):  
  An instance of [`EmptyProblem`](@ref Dionysos.Problem.EmptyProblem) containing the system to abstract and the state region.

- `state_grid` (**required**):  
  The discretization grid for the state space.

- If `state_grid` is not provided, you must set:
    - `h::SVector{N, T}`: Grid spacing vector (used to construct the grid internally).

- `input_grid` (**required**):  
  The discretization grid for the input space.

#### Optional user-tunable fields

- `time_step` (optional, required if the system is continuous):  
  Time step used for discretizing or simulating continuous-time systems.

- `nsystem` (optional, default = `5`):  
  Number of substeps to use when simulating continuous-time systems (e.g., in Runge-Kutta integration).

- `approx_mode` (**required**):  
  The abstraction technique to use. Supported modes:

    - [`USER_DEFINED`](@ref Dionysos.System.DiscreteTimeOverApproximationMap):  
      Use a custom overapproximation function.  
      Set `overapproximation_map::Function`.

    - [`GROWTH`](@ref Dionysos.System.DiscreteTimeGrowthBound):  
      Use growth-bound based overapproximation.  
      Set `jacobian_bound`, optionally `growthbound_map`, `ngrowthbound`.

    - [`LINEARIZED`](@ref Dionysos.System.DiscreteTimeLinearized):  
      Use linearization + Jacobian/Hessian.  
      Set `DF_sys`, `bound_DF`, and `bound_DDF`.

    - [`CENTER_SIMULATION`](@ref Dionysos.System.DiscreteTimeCenteredSimulation):  
      Simulate the center of each cell only.

    - [`RANDOM_SIMULATION`](@ref Dionysos.System.DiscreteTimeRandomSimulation):  
      Sample and simulate random points in each cell.  
      Set `n_samples`.

- `efficient` (optional, default = `true`):  
  Whether to use optimized internal routines based on `approx_mode`.

- `print_level` (optional, default = `1`):  
  Verbosity level:
    - `0`: silent  
    - `1`: standard  
    - `2`: verbose/debug

- `use_periodic_mapping` (optional, default = `false`):  
  If `true`, uses a periodic domain structure when discretizing the state space.

  When enabled, the following fields are required:
    - `periodic_dims::SVector{P, Int}`:  
      Indices of the periodic dimensions.
    - `periodic_periods::SVector{P, T}`:  
      Period length for each periodic dimension.
    - `periodic_start::SVector{P, T}` (optional):  
      Start point of each periodic dimension. Defaults to `0.0` if not provided.

---

### Internally computed fields (after `MOI.optimize!`)

- `abstract_system`:  
  The resulting symbolic abstraction, of type [`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList).

- `discrete_time_system`:  
  Internally constructed discrete-time version of the system used during abstraction.

- `abstraction_construction_time_sec`:  
  Time (in seconds) spent constructing the abstraction.

#### System approximation objects (derived from `approx_mode`)

- `continuous_time_system_approximation`:  
  A [`ContinuousTimeSystemApproximation`](@ref Dionysos.System.ContinuousTimeSystemApproximation), automatically created if the original system is continuous.

- `discrete_time_system_approximation`:  
  A [`DiscreteTimeSystemApproximation`](@ref Dionysos.System.DiscreteTimeSystemApproximation), created in all cases.

### Example

```julia
using Dionysos, JuMP
optimizer = MOI.instantiate(Dionysos.Optim.OptimizerEmptyProblem.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("empty_problem"), my_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), GROWTH)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), my_jacobian_bound)

MOI.optimize!(optimizer)

time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
```
"""
mutable struct OptimizerEmptyProblem{T} <: MOI.AbstractOptimizer
    ## Abstraction Result
    discrete_time_system::Union{
        Nothing,
        MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem,
    }
    abstract_system::Union{Nothing, Dionysos.Symbolic.SymbolicModelList}
    abstraction_construction_time_sec::T

    ## System Approximation
    continuous_time_system_approximation::Union{
        Nothing,
        ST.ContinuousTimeSystemApproximation,
    }
    discrete_time_system_approximation::Union{Nothing, ST.DiscreteTimeSystemApproximation}

    ## User Settings
    empty_problem::Union{Nothing, Dionysos.Problem.EmptyProblem}
    abstraction_region::Any
    incl_mode::MP.INCL_MODE

    ## XMapping & Xset
    XMapping::Union{Nothing, MP.GridMapping}
    state_grid::Union{Nothing, MP.Grid}
    h::Union{Nothing, Any}

    Xset::Union{Nothing, MP.AbstractStateSet}
    Rset::Union{Nothing, MP.AbstractStateSet}

    ## UMapping & Uset
    UMapping::Union{Nothing, MP.GridMapping, MP.ListMapping}

    input_grid::Union{Nothing, MP.Grid}

    Uset::Union{Nothing, MP.AbstractStateSet}

    use_periodic_mapping::Bool
    periodic_dims::Union{Nothing, Any}
    periodic_start::Union{Nothing, Any}
    periodic_periods::Union{Nothing, Any}

    automaton_constructor::Function
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
    threaded::Bool
    efficient::Bool

    print_level::Int
    progress_update_interval::Int
    progress_dt::Float64

    function OptimizerEmptyProblem{T}() where {T}
        optimizer = new{T}(
            nothing,
            nothing,
            0.0,
            nothing,
            nothing,
            nothing,
            nothing,
            MP.INNER,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false, # use periodic domain
            nothing,
            nothing,
            nothing,
            (n, m) -> SY.NewSortedAutomatonList(n, m),
            nothing, #USER_DEFINED
            nothing,
            nothing,
            5, #GROWTH
            nothing,
            nothing,
            nothing, #LINEARIZED
            5, #RANDOM_SIMULATION
            nothing,
            5,
            GROWTH,
            false,
            true,
            1,
            Int(1e5),
            0.2,
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

function reset!(model::OptimizerEmptyProblem)
    model.discrete_time_system = nothing
    model.abstract_system = nothing
    model.abstraction_construction_time_sec = 0.0
    model.continuous_time_system_approximation = nothing
    model.discrete_time_system_approximation = nothing
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

function _pick_state_region(opt::OptimizerEmptyProblem)
    X = opt.abstraction_region
    X === nothing && (X = opt.empty_problem.region)
    X === nothing && (X = opt.empty_problem.system.X)
    return X
end

function build_state_grid(opt::OptimizerEmptyProblem)
    # If user already gave a Grid object, use it directly.
    if opt.state_grid !== nothing
        return opt.state_grid
    end

    # Else build it from h (required).
    if opt.h === nothing
        error("To build the state grid, set either `state_grid` or `h`.")
    end

    if opt.use_periodic_mapping
        _validate_model(opt, [:periodic_dims, :periodic_periods])
        if opt.periodic_start !== nothing
            return MP.get_grid_in_periods(
                opt.periodic_dims,
                opt.periodic_periods,
                opt.periodic_start,
                opt.h,
            )
        else
            return MP.get_grid_in_periods(opt.periodic_dims, opt.periodic_periods, opt.h)
        end
    else
        return MP.GridFree(opt.h)
    end
end

function build_state_mapping(opt::OptimizerEmptyProblem{T}) where {T}
    if opt.XMapping !== nothing
        return opt.XMapping
    end

    grid = build_state_grid(opt)
    X = _pick_state_region(opt)

    # default mapping: explicit enumeration restricted to X
    N = MP.get_dim(grid)
    m = MP.ExplicitGridMapping{N, T}(grid)

    # wrap periodicity if requested
    if opt.use_periodic_mapping
        P = length(opt.periodic_dims)
        start =
            opt.periodic_start === nothing ? SVector{P, T}(ntuple(_->zero(T), P)) :
            opt.periodic_start
        m = MP.PeriodicGridMapping(opt.periodic_dims, opt.periodic_periods, start, m)
    end

    MP.add_set!(m, X, opt.incl_mode)

    return m
end

function build_state_set(opt::OptimizerEmptyProblem)
    if opt.Xset !== nothing
        return opt.Xset
    end
    m = build_state_mapping(opt)
    N = MP.get_dim(m)
    return MP.MappingSet{N}()
end

function build_allowed_state_set(opt::OptimizerEmptyProblem)
    if opt.Rset !== nothing
        return opt.Rset
    end
    return build_state_set(opt)   # Rset === Xset
end

function build_input_mapping(opt::OptimizerEmptyProblem{T}) where {T}
    if opt.UMapping !== nothing
        return opt.UMapping
    end
    M = MP.get_dim(opt.input_grid)
    return MP.ExplicitGridMapping{M, T}(
        opt.input_grid,
        opt.empty_problem.system.U,
        MP.CENTER,
    )
end

function build_input_set(opt::OptimizerEmptyProblem{T}) where {T}
    if opt.Uset !== nothing
        return opt.Uset
    end
    umap = build_input_mapping(opt)
    M = MP.get_dim(umap)
    return MP.MappingSet{M}()
end

_vector_of_tuple(size, value = 0.0) = SVector(ntuple(_ -> value, Val(size)))

function build_noise(optimizer::OptimizerEmptyProblem)
    @warn("Noise is not yet accounted for in system abstraction.")
    concrete_system = optimizer.empty_problem.system
    return _vector_of_tuple(Dionysos.Utils.get_dims(concrete_system.X))
end

function MOI.optimize!(optimizer::OptimizerEmptyProblem)
    t_ref = time()
    # Ensure necessary parameters are set
    _validate_model(optimizer, [:empty_problem])
    @assert optimizer.empty_problem.system !== nothing "System must be set before building overapproximation."

    # Create over-approximation method
    build_system_approximation!(optimizer)

    # Create abstract system
    abstract_system = Dionysos.Symbolic.SymbolicModelList(
        build_state_mapping(optimizer),
        build_input_mapping(optimizer);
        Xset = build_state_set(optimizer),
        Rset = build_allowed_state_set(optimizer),
        Uset = build_input_set(optimizer),
        automaton_constructor = optimizer.automaton_constructor,
    )

    if optimizer.print_level >= 1
        @info("Number of states: $(SY.get_n_state(abstract_system))")
        @info("Number of inputs: $(SY.get_n_input(abstract_system))")
        @info(
            "Number of forward images: $(SY.get_n_input(abstract_system)*SY.get_n_state(abstract_system))"
        )
    end

    # TODO: Consider adding noise handling
    noise = build_noise(optimizer)

    optimizer.print_level >= 1 &&
        println("compute_abstract_system_from_concrete_system!: started")
    if !optimizer.efficient &&
       ST.is_over_approximation(optimizer.discrete_time_system_approximation)
        Dionysos.Symbolic.compute_abstract_system_from_concrete_system!(
            abstract_system,
            ST.get_DiscreteTimeOverApproximationMap(
                optimizer.discrete_time_system_approximation,
            );
            verbose = optimizer.print_level >= 2,
            update_interval = optimizer.progress_update_interval,
            progress_dt = optimizer.progress_dt,
            threaded = optimizer.threaded,
        )
    else
        Dionysos.Symbolic.compute_abstract_system_from_concrete_system!(
            abstract_system,
            optimizer.discrete_time_system_approximation;
            verbose = optimizer.print_level >= 2,
            update_interval = optimizer.progress_update_interval,
            progress_dt = optimizer.progress_dt,
            threaded = optimizer.threaded,
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
