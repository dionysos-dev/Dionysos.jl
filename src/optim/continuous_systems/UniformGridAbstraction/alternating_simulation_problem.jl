@enum ApproxMode USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION

"""
    OptimizerAlternatingSimulationProblem{T} <: Dionysos.Optim.AbstractDionysosOptimizer

A solver responsible for constructing a **symbolic abstraction of the system dynamics**,
independently of any control specification.

This optimizer wraps everything needed to solve an [`AlternatingSimulationProblem`](@ref Dionysos.Problem.AlternatingSimulationProblem),
which is used to generate a symbolic model (abstraction) of either a continuous- or discrete-time system.

The optimizer supports several abstraction modes, optional implicit mappings/state sets,
periodic mappings, multithreaded computation, and distributed partition-based abstraction.

---

### Purpose

This optimizer builds a symbolic abstraction by discretizing the state and input spaces,
constructing a state/input mapping, defining a source-state domain, and computing the
abstract transition relation from a system approximation.

The abstraction method is selected via the `approx_mode` field. Depending on the chosen mode,
the optimizer constructs either an over-approximation, an under-approximation, or a simulation-based abstraction.

The resulting abstract model can then be reused by higher-level control solvers
(safety, reachability, co-safe LTL, etc.).

---

### Abstraction semantics

The constructed symbolic model distinguishes between:

- `XMapping`: global mapping from concrete states to abstract states.
- `Xset`: set of **source states** that are enumerated when building transitions.
- `Rset`: set of **retained / allowed states** that may appear as transition targets.
- `UMapping`: mapping from concrete inputs to abstract inputs.
- `Uset`: admissible abstract-input set.

In the standard non-distributed setting, `Xset` usually coincides with `Rset`.
In distributed mode, the abstraction may be computed on local partitions of `Xset`,
while retaining a common global `Rset`.

---

### Parameters

#### Mandatory fields set by the user

- `alternating_simulation_problem` (**required**):  
  An instance of [`AlternatingSimulationProblem`](@ref Dionysos.Problem.AlternatingSimulationProblem) containing the system to abstract and the target abstraction region.

- State discretization (**required**):
    - either `state_grid`
    - or `h`, from which the state grid is built internally.

- `input_grid` (**required**, unless `UMapping` is set directly):  
  The discretization grid for the input space.

#### Optional mapping / set fields

- `abstraction_region` (optional):  
  Concrete region used to define the abstraction domain.  
  If not provided, the optimizer uses:
  1. `alternating_simulation_problem.region`, if available,
  2. otherwise `alternating_simulation_problem.system.X`.

- `incl_mode` (optional, default = `MP.INNER`):  
  Inclusion mode used when constructing mappings or state sets.

- `XMapping` (optional):  
  Pre-built abstract-state mapping. If not provided, it is constructed from the state grid.

- `Xset` (optional):  
  Source abstract-state set. If not provided, it is built from `abstraction_region`.

- `Rset` (optional):  
  Retained / allowed target-state set. If not provided, it defaults to `copy(Xset)`.

- `UMapping` (optional):  
  Pre-built abstract-input mapping.

- `Uset` (optional):  
  Abstract-input set. If not provided, a default admissible set is built.

#### Grid / implicit representation settings

- `state_grid` (optional):  
  Explicit state grid used for discretization.

- `h` (optional):  
  Grid spacing vector used to construct the state grid if `state_grid` is not provided.

- `use_implicit_mapping` (optional, default = `false`):  
  If `true`, constructs an implicit state mapping instead of an explicit one.

- `mapping_region` (required if `use_implicit_mapping = true`):  
  Hyper-rectangle used as ambient region for the implicit mapping.  
  It must enclose `abstraction_region`.

- `use_implicit_stateset` (optional, default = `false`):  
  If `true`, builds `Xset` as an implicit state set rather than an explicit set of indices.

#### Periodic mapping settings

- `use_periodic_mapping` (optional, default = `false`):  
  If `true`, wraps the state mapping as a periodic mapping.

When enabled, the following fields are required:

- `periodic_dims`: indices of periodic dimensions.
- `periodic_periods`: period length for each periodic dimension.
- `periodic_start` (optional): start point for each periodic dimension. Defaults to zero if not provided.

#### System approximation settings

- `approx_mode` (optional, default = `GROWTH`):  
  Abstraction technique used to build the system approximation. Supported modes:

    - [`USER_DEFINED`](@ref Dionysos.System.DiscreteTimeOverApproximationMap):  
      Use a custom over-approximation map.  
      Set `overapproximation_map`.

    - [`GROWTH`](@ref Dionysos.System.DiscreteTimeGrowthBound):  
      Use growth-bound based approximation.  
      Set `jacobian_bound`, or directly `growthbound_map`.  
      `ngrowthbound` controls the internal growth-bound discretization parameter.

    - [`LINEARIZED`](@ref Dionysos.System.DiscreteTimeLinearized):  
      Use linearization and derivative bounds.  
      Set `DF_sys`, `bound_DF`, and `bound_DDF`.

    - [`CENTER_SIMULATION`](@ref Dionysos.System.DiscreteTimeCenteredSimulation):  
      Simulate the center of each abstract cell only.

    - [`RANDOM_SIMULATION`](@ref Dionysos.System.DiscreteTimeRandomSimulation):  
      Simulate randomly sampled points in each abstract cell.  
      Set `n_samples`.

- `efficient` (optional, default = `true`):
  Deprecated and ignored — the abstraction kernel hoists approximation-specific
  per-input work uniformly (see `Dionysos.System.input_cache` / `reach_set`).

#### Continuous-time settings

- `time_step` (required for continuous-time systems):  
  Sampling step used to discretize a continuous-time approximation.

- `nsystem` (optional, default = `5`):  
  Number of internal substeps used during continuous-time simulation/discretization routines.

#### Execution settings

- `execution_backend` (optional, default = `SY.SequentialBackend()`):  
  Execution backend used to compute the transition relation.

Supported execution backends:

- `SY.SequentialBackend()`  
  Sequential computation on the current Julia process.

- `SY.ThreadedBackend(progress_dt)`  
  Multithreaded computation on the current Julia process.

- `SY.JuliaDistributedBackend(procs, nparts, partition_strategy, threaded_per_worker)`
  Distributed computation over Julia worker processes.

- `SY.SlurmArrayBackend(nchunks, chunk_id, outdir, partition_strategy, write_only)`  
  SLURM-array execution. Each array task computes one chunk and writes it to disk.

#### Logging / progress settings

- `print_level` (optional, default = `1`):  
  Verbosity level:
    - `0`: silent
    - `1`: summary information
    - `2`: detailed progress output

- `progress_update_interval` (optional, default = `Int(1e5)`):  
  Number of source-state/input pairs processed between progress updates.

- `progress_dt` (optional, default = `0.2`):  
  Minimum wall-clock time between progress refreshes in threaded mode.

---

### Internally computed fields (after `MOI.optimize!`)

- `abstract_system`:  
  The resulting symbolic abstraction, of type [`SymbolicModelList`](@ref SY.SymbolicModelList).

- `discrete_time_system`:  
  Internally constructed discrete-time system used to generate the abstraction.

- `abstraction_construction_time_sec`:  
  Total abstraction-construction time in seconds.

#### Approximation objects derived from `approx_mode`

- `continuous_time_system_approximation`:  
  A [`ContinuousTimeSystemApproximation`](@ref Dionysos.System.ContinuousTimeSystemApproximation), built automatically if the original system is continuous.

- `discrete_time_system_approximation`:  
  A [`DiscreteTimeSystemApproximation`](@ref Dionysos.System.DiscreteTimeSystemApproximation), used by the abstraction kernel in all cases.

---

### Typical workflow

1. Define the concrete empty problem and discretization parameters.
2. Build or infer the state/input mappings and state sets.
3. Construct the system approximation from `approx_mode`.
4. Compute the abstract transition relation, either:
    - sequentially,
    - multithreaded,
    - or distributed across source-state partitions.
5. Store the resulting symbolic abstraction in `abstract_system`.

---

### Example

```julia
using Dionysos, JuMP
optimizer = MOI.instantiate(Dionysos.Optim.OptimizerAlternatingSimulationProblem.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("alternating_simulation_problem"), my_problem)
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
mutable struct OptimizerAlternatingSimulationProblem{T} <: OP.AbstractDionysosOptimizer
    ## Abstraction Result
    discrete_time_system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstraction_construction_time_sec::T

    ## System Approximation
    continuous_time_system_approximation::Union{
        Nothing,
        ST.ContinuousTimeSystemApproximation,
    }
    discrete_time_system_approximation::Union{Nothing, ST.DiscreteTimeSystemApproximation}

    ## User Settings
    alternating_simulation_problem::Union{
        Nothing,
        Dionysos.Problem.AlternatingSimulationProblem,
    }
    abstraction_region::Any
    incl_mode::MP.INCL_MODE

    ## XMapping & Xset
    XMapping::Union{Nothing, MP.GridMapping}
    state_grid::Union{Nothing, MP.Grid}
    h::Union{Nothing, Any}
    use_implicit_mapping::Bool
    mapping_region::Union{Nothing, UT.Box}

    Xset::Union{Nothing, MP.AbstractStateSet}
    Rset::Union{Nothing, MP.AbstractStateSet}
    use_implicit_stateset::Bool

    state_filter::Union{Nothing, Function}
    state_input_filter::Union{Nothing, Function}

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

    ### Execution settings
    execution_backend::SY.AbstractExecutionBackend
    approx_mode::ApproxMode
    efficient::Bool

    ### Metadata settings
    transition_metadata::SY.AbstractTransitionMetadata

    print_level::Int
    progress_update_interval::Int
    progress_dt::Float64

    function OptimizerAlternatingSimulationProblem{T}() where {T}
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
            false, # use implicit mapping
            nothing,
            nothing,
            nothing,
            false, # use implicit stateset
            nothing, # state_filter 
            nothing, # state_input_filter
            nothing,
            nothing,
            nothing,
            false, # use periodic domain
            nothing,
            nothing,
            nothing,
            (n, m) -> SY.SortedAutomatonList(n, m),
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
            SY.SequentialBackend(), # execution_backend
            GROWTH,                 # approx
            true,                   # efficient
            SY.NoTransitionMetadata(),
            1,
            Int(1e5),
            0.2,
        )
        return optimizer
    end
end

OptimizerAlternatingSimulationProblem() = OptimizerAlternatingSimulationProblem{Float64}()

MOI.is_empty(optimizer::OptimizerAlternatingSimulationProblem) =
    optimizer.alternating_simulation_problem === nothing

function MOI.get(model::OptimizerAlternatingSimulationProblem, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec
end

function reset!(model::OptimizerAlternatingSimulationProblem)
    model.discrete_time_system = nothing
    model.abstract_system = nothing
    model.abstraction_construction_time_sec = 0.0
    model.continuous_time_system_approximation = nothing
    model.discrete_time_system_approximation = nothing
    return model
end

function _validate_model(
    model::OptimizerAlternatingSimulationProblem,
    required_fields::Vector{Symbol},
)
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set the `$(field)`. Missing required field in OptimizerAlternatingSimulationProblem.",
            )
        end
    end
end

function build_continuous_approximation(
    optimizer::OptimizerAlternatingSimulationProblem,
    system,
)
    mode = optimizer.approx_mode
    if mode == USER_DEFINED
        _validate_model(optimizer, [:overapproximation_map])
        return ST.ContinuousTimeOverApproximationMap(
            system,
            optimizer.overapproximation_map,
        )
    elseif mode == GROWTH
        if optimizer.growthbound_map !== nothing
            return ST.ContinuousTimeGrowthBound(system, optimizer.growthbound_map)
        else
            return ST.ContinuousTimeGrowthBound(
                system;
                jacobian_bound = optimizer.jacobian_bound,
                ngrowthbound = optimizer.ngrowthbound,
            )
        end
    elseif mode == LINEARIZED
        _validate_model(optimizer, [:DF_sys, :bound_DF, :bound_DDF])
        return ST.ContinuousTimeLinearized(
            system,
            optimizer.DF_sys,
            optimizer.bound_DF,
            optimizer.bound_DDF;
            num_substeps = optimizer.nsystem,
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

function build_discrete_approximation(
    optimizer::OptimizerAlternatingSimulationProblem,
    system,
)
    mode = optimizer.approx_mode
    if mode == USER_DEFINED
        _validate_model(optimizer, [:overapproximation_map])
        return ST.DiscreteTimeOverApproximationMap(system, optimizer.overapproximation_map)
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

function build_system_approximation!(optimizer::OptimizerAlternatingSimulationProblem)
    _validate_model(optimizer, [:alternating_simulation_problem])
    @assert optimizer.alternating_simulation_problem.system !== nothing "System must be set before building overapproximation."

    system = optimizer.alternating_simulation_problem.system
    if isa(system, MS.ConstrainedBlackBoxControlContinuousSystem)
        _validate_model(optimizer, [:time_step])  # Ensure time step is provided
        optimizer.continuous_time_system_approximation =
            build_continuous_approximation(optimizer, system)
        optimizer.discrete_time_system_approximation = ST.discretize(
            optimizer.continuous_time_system_approximation,
            optimizer.time_step;
            num_substeps = optimizer.nsystem,
        )
    elseif isa(system, MS.ConstrainedBlackBoxControlDiscreteSystem)
        optimizer.discrete_time_system_approximation =
            build_discrete_approximation(optimizer, system)
    else
        error("Unknown system type: $(typeof(system))")
    end
    return optimizer.discrete_time_system =
        ST.get_system(optimizer.discrete_time_system_approximation)
end

function _validate_periodic_data(opt)
    opt.use_periodic_mapping || return

    pd = opt.periodic_dims
    pp = opt.periodic_periods
    ps = opt.periodic_start

    P = length(pd)
    length(pp) == P || error("periodic_periods must have length $P, got $(length(pp))")
    length(ps) == P || error("periodic_start must have length $P, got $(length(ps))")

    N = LazySets.dim(opt.alternating_simulation_problem.system.X)
    if P > 0
        all(1 .<= pd .<= N) || error("periodic_dims must be in 1:$N, got $pd")
    end
end

function build_abstraction_region(opt::OptimizerAlternatingSimulationProblem)
    X = opt.abstraction_region
    X === nothing && (X = opt.alternating_simulation_problem.region)
    X === nothing && (X = opt.alternating_simulation_problem.system.X)
    return X
end

function build_state_grid(opt::OptimizerAlternatingSimulationProblem)
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
        _validate_periodic_data(opt)
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

function _validate_mapping_encloses_abstraction_region(opt)
    opt.use_implicit_mapping || return

    _validate_model(opt, [:mapping_region, :abstraction_region])

    maprect = opt.mapping_region
    absbox = UT._outer_box(opt.abstraction_region)

    if !Base.issubset(absbox, maprect)
        error(
            "mapping_region must fully enclose abstraction_region.\n" *
            "mapping_region = $maprect\n" *
            "outer box of abstraction_region = $absbox",
        )
    end
end

function build_state_mapping(opt::OptimizerAlternatingSimulationProblem{T}) where {T}
    if opt.XMapping !== nothing
        return opt.XMapping
    end

    grid = build_state_grid(opt)
    if opt.use_implicit_mapping
        _validate_model(opt, [:mapping_region, :abstraction_region])
        _validate_mapping_encloses_abstraction_region(opt)
        maprect = opt.mapping_region
        m = MP.ImplicitGridMapping(grid, maprect; incl_mode = opt.incl_mode)
    else
        N = MP.get_dim(grid)
        m = MP.ExplicitGridMapping{N, T}(grid)
    end

    if opt.use_periodic_mapping
        P = length(opt.periodic_dims)
        start =
            opt.periodic_start === nothing ? SVector{P, T}(ntuple(_ -> zero(T), P)) :
            opt.periodic_start
        m = MP.PeriodicGridMapping(opt.periodic_dims, opt.periodic_periods, start, m)
    end

    if !opt.use_implicit_mapping
        MP.add_set!(m, opt.abstraction_region, opt.incl_mode)
    end

    return m
end

function build_state_set(opt::OptimizerAlternatingSimulationProblem)
    if opt.Xset !== nothing
        return opt.Xset
    end

    XMapping = opt.XMapping
    if opt.use_implicit_stateset
        S = MP.ImplicitStateSet(XMapping, opt.abstraction_region, opt.incl_mode)
    else
        S = MP.ExplicitIdSet{MP.get_dim(XMapping)}()
        MP.add_set!(S, XMapping, opt.abstraction_region, opt.incl_mode)
    end
    return S
end

function build_allowed_state_set(opt::OptimizerAlternatingSimulationProblem)
    return opt.Rset === nothing ? copy(opt.Xset) : opt.Rset
end

function build_input_mapping(opt::OptimizerAlternatingSimulationProblem{T}) where {T}
    if opt.UMapping !== nothing
        return opt.UMapping
    end
    M = MP.get_dim(opt.input_grid)
    return MP.ExplicitGridMapping{M, T}(
        opt.input_grid,
        opt.alternating_simulation_problem.system.U,
        MP.CENTER,
    )
end

function build_input_set(opt::OptimizerAlternatingSimulationProblem{T}) where {T}
    if opt.Uset !== nothing
        return opt.Uset
    end
    umap = build_input_mapping(opt)
    M = MP.get_dim(umap)
    return MP.MappingSet{M}()
end

_vector_of_tuple(size, value = 0.0) = SVector(ntuple(_ -> value, Val(size)))

function build_noise(optimizer::OptimizerAlternatingSimulationProblem)
    @warn("Noise is not yet accounted for in system abstraction.")
    concrete_system = optimizer.alternating_simulation_problem.system
    return _vector_of_tuple(LazySets.dim(concrete_system.X))
end

function build_empty_abstraction!(optimizer::OptimizerAlternatingSimulationProblem)
    _validate_model(optimizer, [:alternating_simulation_problem])

    build_system_approximation!(optimizer)

    optimizer.abstraction_region = build_abstraction_region(optimizer)
    optimizer.XMapping = build_state_mapping(optimizer)
    optimizer.UMapping = build_input_mapping(optimizer)
    optimizer.Xset = build_state_set(optimizer)
    optimizer.Rset = build_allowed_state_set(optimizer)
    optimizer.Uset = build_input_set(optimizer)

    optimizer.abstract_system = SY.SymbolicModelList(
        optimizer.XMapping,
        optimizer.UMapping;
        Xset = optimizer.Xset,
        Rset = optimizer.Rset,
        Uset = optimizer.Uset,
        automaton_constructor = optimizer.automaton_constructor,
        metadata = optimizer.transition_metadata,
    )

    return optimizer.abstract_system
end

function MOI.optimize!(optimizer::OptimizerAlternatingSimulationProblem)
    t_ref = time()

    _validate_model(optimizer, [:alternating_simulation_problem])
    @assert optimizer.alternating_simulation_problem.system !== nothing

    abstract_system = build_empty_abstraction!(optimizer)

    if optimizer.print_level >= 1
        @info("Number of states: $(SY.get_n_state(abstract_system))")
        @info("Number of inputs: $(SY.get_n_input(abstract_system))")
        @info(
            "Number of forward images: $(SY.get_n_input(abstract_system) * SY.get_n_state(abstract_system))",
        )
    end

    build_noise(optimizer)

    optimizer.print_level >= 1 &&
        println("compute_abstract_system_from_concrete_system!: started")

    SY.compute_abstract_system_from_concrete_system!(
        abstract_system,
        optimizer.discrete_time_system_approximation;
        execution_backend = optimizer.execution_backend,
        print_level = optimizer.print_level,
        update_interval = optimizer.progress_update_interval,
        state_filter = optimizer.state_filter,
        state_input_filter = optimizer.state_input_filter,
    )

    optimizer.print_level >= 1 && println(
        "compute_abstract_system_from_concrete_system! terminated with success: ",
        "$(HybridSystems.ntransitions(SY.get_automaton(abstract_system))) transitions created",
    )

    optimizer.abstract_system = abstract_system
    optimizer.abstraction_construction_time_sec = time() - t_ref
    return
end
