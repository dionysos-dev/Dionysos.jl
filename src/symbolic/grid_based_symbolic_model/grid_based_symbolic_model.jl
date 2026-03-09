"""
    GridBasedSymbolicModel{N,M} <: SymbolicModel{N,M}

Intermediate abstract type for symbolic models that rely on a grid-based or
mapping-based discretization.

Semantics:
- state mapping: global abstract-state numbering / coordinate map
- input mapping: global abstract-input numbering / coordinate map
- source domain (`Xset`): states enumerated as sources
- retained domain (`Rset`): states allowed as targets
"""
abstract type GridBasedSymbolicModel{N,M} <: SymbolicModel{N,M} end

# ------------------------------------------------------------------
# Required interface for concrete subtypes
# ------------------------------------------------------------------

get_source_domain(symmodel::GridBasedSymbolicModel) =
    error("get_source_domain not implemented for $(typeof(symmodel))")

get_retained_domain(symmodel::GridBasedSymbolicModel) =
    error("get_retained_domain not implemented for $(typeof(symmodel))")

# ------------------------------------------------------------------
# Default generic behavior
# ------------------------------------------------------------------

enum_source_states(symmodel::GridBasedSymbolicModel) =
    MP.enum_states(get_source_domain(symmodel), get_state_mapping(symmodel))

# enum_inputs(symmodel::GridBasedSymbolicModel) =
#     MP.enum_inputs(get_input_mapping(symmodel))

is_allowed_state(symmodel::GridBasedSymbolicModel, q) =
    q !== nothing &&
    MP.contains_state(get_retained_domain(symmodel), get_state_mapping(symmodel), q)

get_n_state(symmodel::GridBasedSymbolicModel) =
    length(collect(enum_source_states(symmodel)))

get_n_input(symmodel::GridBasedSymbolicModel) =
    length(collect(enum_inputs(symmodel)))

get_state_domain(symmodel::GridBasedSymbolicModel) = get_source_domain(symmodel)
enum_states(symmodel::GridBasedSymbolicModel) = enum_source_states(symmodel)

# ------------------------------------------------
# Transition collection
# ------------------------------------------------

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::Union{
        ST.DiscreteTimeSystemOverApproximation,
        ST.DiscreteTimeGrowthBound,
        ST.DiscreteTimeLinearized,
        ST.DiscreteTimeSystemUnderApproximation,
        ST.DiscreteTimeCenteredSimulation,
    };
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    trans = collect_abstract_transitions(
        abstract_system,
        concrete_system_approx;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
    isempty(trans) || add_transitions!(abstract_system, trans)
    return
end

function collect_abstract_transitions(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    trans = Tuple{Int,Int,Int}[]
    collect_abstract_transitions!(
        trans,
        symmodel,
        concrete_system_approx;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
    return trans
end

# ------------------------------------------------
# Low-level transition helpers
# ------------------------------------------------

function compute_abstract_transitions_from_rectangle!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::UT.HyperRectangle,
    abstract_state::Int,
    abstract_input::Int,
    translist::Vector{Tuple{Int,Int,Int}},
)
    targets, allin = get_states_from_set_strict(symmodel, reachable_set, MP.OUTER)
    allin || return false

    for target in targets
        is_allowed_state(symmodel, target) || return false
    end

    for target in targets
        push!(translist, (target, abstract_state, abstract_input))
    end

    return true
end

function compute_abstract_transitions_from_points!(
    symmodel::GridBasedSymbolicModel,
    reachable_points,
    abstract_state::Int,
    abstract_input::Int,
    translist::Vector{Tuple{Int,Int,Int}},
)
    start_len = length(translist)

    for y in reachable_points
        target = get_abstract_state(symmodel, y)
        if target === nothing || !is_allowed_state(symmodel, target)
            resize!(translist, start_len)
            return false
        end
        push!(translist, (target, abstract_state, abstract_input))
    end

    unique!(view(translist, start_len+1:length(translist)))
    unique!(translist)  # simple first version
    return true
end


# ------------------------------------------------
# Generic helpers for sequential/threaded execution
# ------------------------------------------------

# Unified Dispatcher
function _collect_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    if !threaded || Threads.nthreads() == 1
        return _collect_transitions_sequential!(
            out,
            symmodel,
            workfun!;
            verbose = verbose,
            update_interval = update_interval,
        )
    else
        return _collect_transitions_threaded!(
            out,
            symmodel,
            workfun!;
            verbose = verbose,
            update_interval = update_interval,
            progress_dt = progress_dt,
        )
    end
end

# Sequential double loop
function _collect_transitions_sequential!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
)
    inputs = collect(enum_inputs(symmodel))
    states = collect(enum_source_states(symmodel))

    total_work = length(inputs) * length(states)
    total_updates = max(div(total_work, max(1, update_interval)), 1)
    progress = verbose ? ProgressMeter.Progress(total_updates) : nothing

    localbuf = Tuple{Int,Int,Int}[]
    count = 0

    for abstract_input in inputs
        for abstract_state in states
            empty!(localbuf)
            workfun!(localbuf, abstract_state, abstract_input)
            append!(out, localbuf)

            count += 1
            if verbose && count % update_interval == 0
                ProgressMeter.next!(progress)
            end
        end
    end

    verbose && ProgressMeter.finish!(progress)
    return out
end

# Threaded double loop
function _collect_transitions_threaded!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
)
    inputs = collect(enum_inputs(symmodel))
    states = collect(enum_source_states(symmodel))

    total_work = length(inputs) * length(states)
    nthreads = Threads.nthreads()

    transitions_by_thread = [Vector{Tuple{Int,Int,Int}}() for _ in 1:nthreads]
    local_done = fill(0, nthreads)

    progress_dt_ns = Int(round(progress_dt * 1e9))
    prog = verbose ? ProgressMeter.Progress(total_work) : nothing
    global_done = Threads.Atomic{Int}(0)
    last_t = time_ns()

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        abstract_input = @inbounds inputs[input_idx]
        abstract_state = @inbounds states[state_idx]

        prev_length = length(local_transitions)
        workfun!(local_transitions, abstract_state, abstract_input)

        local_done[tid] += 1
        if local_done[tid] >= update_interval
            Threads.atomic_add!(global_done, local_done[tid])
            local_done[tid] = 0
        end

        if verbose && tid == 1
            t = time_ns()
            if t - last_t >= progress_dt_ns
                ProgressMeter.update!(prog, global_done[] + local_done[1])
                last_t = t
            end
        end
    end

    for tid in 1:nthreads
        if local_done[tid] > 0
            Threads.atomic_add!(global_done, local_done[tid])
        end
    end

    if verbose
        ProgressMeter.update!(prog, global_done[])
        ProgressMeter.finish!(prog)
    end

    for local_transitions in transitions_by_thread
        isempty(local_transitions) || append!(out, local_transitions)
    end

    return out
end

# ------------------------------------------------
# Approximation Kernels
# ------------------------------------------------

# Kernel that computes transitions from a reachable set over-approximation
function collect_abstract_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemOverApproximation;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    compute_reachable_set = ST.get_over_approximation_map(concrete_system_approx)

    workfun! = function (
        transbuf::Vector{Tuple{Int,Int,Int}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = get_concrete_input(symmodel, abstract_input)
        concrete_elem = get_concrete_elem(symmodel, abstract_state)
        reachable_set = compute_reachable_set(concrete_elem, concrete_input)

        localbuf = Tuple{Int,Int,Int}[]
        allin = compute_abstract_transitions_from_rectangle!(
            symmodel,
            reachable_set,
            abstract_state,
            abstract_input,
            localbuf,
        )
        allin && append!(transbuf, localbuf)
        return nothing
    end

    return _collect_transitions!(
        out,
        symmodel,
        workfun!;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
end

# Kernel that computes transitions from a reachable set over-approximation
# using growth bound
function collect_abstract_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeGrowthBound;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    growthbound_map = concrete_system_approx.growthbound_map
    system_map = ST.get_system_map(concrete_system_approx)

    XMapping = get_state_mapping(symmodel)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0

    inputs = collect(enum_inputs(symmodel))
    input_data = Dict{Int,Tuple{Any,Any}}()
    for abstract_input in inputs
        concrete_input = get_concrete_input(symmodel, abstract_input)
        Fr = growthbound_map(r, concrete_input)
        input_data[abstract_input] = (concrete_input, Fr)
    end

    workfun! = function (
        transbuf::Vector{Tuple{Int,Int,Int}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input, Fr = input_data[abstract_input]
        concrete_state = get_concrete_state(symmodel, abstract_state)
        Fx = system_map(concrete_state, concrete_input)
        reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)

        localbuf = Tuple{Int,Int,Int}[]
        allin = compute_abstract_transitions_from_rectangle!(
            symmodel,
            reachable_set,
            abstract_state,
            abstract_input,
            localbuf,
        )
        allin && append!(transbuf, localbuf)
        return nothing
    end

    return _collect_transitions!(
        out,
        symmodel,
        workfun!;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
end

# Kernel that computes transitions from a reachable set over-approximation
# using linearization
function collect_abstract_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeLinearized;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    XMapping = get_state_mapping(symmodel)
    N = MP.get_dim(XMapping)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0
    _H_ = SMatrix{N,N}(LinearAlgebra.I) .* r
    _ONE_ = ones(SVector{N})
    e = LinearAlgebra.norm(r, Inf)

    error_map = concrete_system_approx.error_map
    linsys_map = concrete_system_approx.linsys_map

    inputs = collect(enum_inputs(symmodel))
    input_data = Dict{Int,Tuple{Any,Any,Any}}()
    for abstract_input in inputs
        concrete_input = get_concrete_input(symmodel, abstract_input)
        Fe = error_map(e, concrete_input)
        Fr = r .+ Fe
        input_data[abstract_input] = (concrete_input, Fe, Fr)
    end

    workfun! = function (
        transbuf::Vector{Tuple{Int,Int,Int}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input, Fe, Fr = input_data[abstract_input]
        concrete_state = get_concrete_state(symmodel, abstract_state)

        Fx, DFx = linsys_map(concrete_state, _H_, concrete_input)

        A = inv(DFx)
        b = abs.(A) * Fr .+ 1.0
        HP = UT.CenteredPolyhedron(A, b)

        rad = abs.(DFx) * _ONE_ .+ Fe
        reachable_set = UT.HyperRectangle(Fx - rad, Fx + rad)

        localbuf = Tuple{Int,Int,Int}[]
        allin = compute_abstract_transitions_from_rectangle!(
            symmodel,
            reachable_set,
            abstract_state,
            abstract_input,
            localbuf,
        )
        allin && append!(transbuf, localbuf)
        return nothing
    end

    return _collect_transitions!(
        out,
        symmodel,
        workfun!;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
end

# Kernel that computes transitions from a reachable set under-approximation
function collect_abstract_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemUnderApproximation;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    under_approximation_map = ST.get_under_approximation_map(concrete_system_approx)

    workfun! = function (
        transbuf::Vector{Tuple{Int,Int,Int}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = get_concrete_input(symmodel, abstract_input)
        concrete_elem = get_concrete_elem(symmodel, abstract_state)
        reachable_points = under_approximation_map(concrete_elem, concrete_input)

        localbuf = Tuple{Int,Int,Int}[]
        allin = compute_abstract_transitions_from_points!(
            symmodel,
            reachable_points,
            abstract_state,
            abstract_input,
            localbuf,
        )
        allin && append!(transbuf, localbuf)
        return nothing
    end

    return _collect_transitions!(
        out,
        symmodel,
        workfun!;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
end

# Kernel that computes transitions from a reachable set under-approximation
# using center simulation
function collect_abstract_transitions!(
    out::Vector{Tuple{Int,Int,Int}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeCenteredSimulation;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    system_map = ST.get_system_map(concrete_system_approx)

    workfun! = function (
        transbuf::Vector{Tuple{Int,Int,Int}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = get_concrete_input(symmodel, abstract_input)
        concrete_state = get_concrete_state(symmodel, abstract_state)
        y = system_map(concrete_state, concrete_input)
        target = get_abstract_state(symmodel, y)

        if target !== nothing && is_allowed_state(symmodel, target)
            push!(transbuf, (target, abstract_state, abstract_input))
        end
        return nothing
    end

    return _collect_transitions!(
        out,
        symmodel,
        workfun!;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )
end

# ------------------------------------------------
# Support for Muliprocessing
# ------------------------------------------------
import Distributed

# ------------------------------------------------
# Local partition model wrapper
# ------------------------------------------------

"""
    LocalGridBasedSymbolicModel

Wrapper around a global symbolic model that overrides only the source-domain
enumeration while keeping the same global state/input mappings and retained set.
"""
struct LocalGridBasedSymbolicModel{N,M,SM,XS} <: GridBasedSymbolicModel{N,M}
    parent::SM
    Xset_local::XS
end

LocalGridBasedSymbolicModel(symmodel::GridBasedSymbolicModel{N,M}, Xset_local) where {N,M} =
    LocalGridBasedSymbolicModel{N,M,typeof(symmodel),typeof(Xset_local)}(symmodel, Xset_local)

get_state_mapping(sym::LocalGridBasedSymbolicModel) = get_state_mapping(sym.parent)
get_input_mapping(sym::LocalGridBasedSymbolicModel) = get_input_mapping(sym.parent)
get_source_domain(sym::LocalGridBasedSymbolicModel) = sym.Xset_local
get_retained_domain(sym::LocalGridBasedSymbolicModel) = get_retained_domain(sym.parent)

# forward everything else to parent
get_concrete_state(sym::LocalGridBasedSymbolicModel, q) = get_concrete_state(sym.parent, q)
get_concrete_elem(sym::LocalGridBasedSymbolicModel, q) = get_concrete_elem(sym.parent, q)
get_concrete_input(sym::LocalGridBasedSymbolicModel, u) = get_concrete_input(sym.parent, u)
get_abstract_state(sym::LocalGridBasedSymbolicModel, x) = get_abstract_state(sym.parent, x)
add_transitions!(sym::LocalGridBasedSymbolicModel, trans) = add_transitions!(sym.parent, trans)

# ------------------------------------------------------------------
# Partitioning
# ------------------------------------------------------------------

"""
    partition_source_states(symmodel, nparts; strategy=:roundrobin)

Partition the source states of `symmodel` into `nparts` subsets.

Currently supported strategies:
- `:roundrobin`
- `:contiguous`
"""
function partition_source_states(
    symmodel::GridBasedSymbolicModel,
    nparts::Int;
    strategy::Symbol = :roundrobin,
)
    nparts >= 1 || error("nparts must be >= 1")

    states = collect(enum_source_states(symmodel))
    isempty(states) && return [ExplicitStateSubset(Int[]) for _ in 1:nparts]

    parts = [Int[] for _ in 1:nparts]

    if strategy == :roundrobin
        for (k, q) in enumerate(states)
            push!(parts[mod1(k, nparts)], q)
        end
    elseif strategy == :contiguous
        n = length(states)
        for i in 1:nparts
            a = fld((i - 1) * n, nparts) + 1
            b = fld(i * n, nparts)
            if a <= b
                append!(parts[i], @view states[a:b])
            end
        end
    else
        error("Unknown partition strategy: $strategy")
    end

    return [MP.ExplicitIdSet(p) for p in parts]
end

# ------------------------------------------------------------------
# Local worker kernel
# ------------------------------------------------------------------

"""
    _run_local_partition(symmodel, concrete_system_approx; kwargs...)

Compute transitions for one local partition only.
"""
function _run_local_partition(
    local_symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    transitions = collect_abstract_transitions(
        local_symmodel,
        concrete_system_approx;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )

    n_source_states = length(collect(enum_source_states(local_symmodel)))

    return DistributedAbstractionResult(
        transitions,
        n_source_states,
        length(transitions),
    )
end

# ------------------------------------------------------------------
# Public distributed collector
# ------------------------------------------------------------------

"""
    collect_abstract_transitions_distributed(symmodel, concrete_system_approx; kwargs...)

Compute abstract transitions by partitioning the source states across multiple
Julia worker processes.

This function returns the transition tuples but does not mutate `symmodel`.
"""
function collect_abstract_transitions_distributed(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    procs = Distributed.workers(),
    nparts::Int = length(procs),
    partition_strategy::Symbol = :roundrobin,
    threaded_per_worker::Bool = false,
    verbose::Bool = true,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
)
    isempty(procs) && error("No distributed workers available. Call Distributed.addprocs(...) first.")
    nparts >= 1 || error("nparts must be >= 1")

    parts = partition_source_states(symmodel, nparts; strategy = partition_strategy)
    jobs = [(LocalGridBasedSymbolicModel(symmodel, Xset_local), concrete_system_approx) for Xset_local in parts]

    # Keep worker-side verbosity off; otherwise outputs interleave badly.
    results = Distributed.pmap(jobs) do job
        local_symmodel, local_approx = job
        _run_local_partition(
            local_symmodel,
            local_approx;
            verbose = false,
            update_interval = update_interval,
            progress_dt = progress_dt,
            threaded = threaded_per_worker,
        )
    end

    transitions = Tuple{Int,Int,Int}[]
    total_sources = 0
    total_transitions = 0

    for res in results
        append!(transitions, res.transitions)
        total_sources += res.n_source_states
        total_transitions += res.n_transitions
    end

    verbose && @info(
        "Distributed abstraction finished",
        nparts = nparts,
        total_source_states = total_sources,
        total_transitions = total_transitions,
    )

    return transitions
end

# ------------------------------------------------------------------
# Public distributed mutating API
# ------------------------------------------------------------------

"""
    compute_abstract_system_distributed!(symmodel, concrete_system_approx; kwargs...)

Distributed version of `compute_abstract_system_from_concrete_system!`.

It computes transitions in parallel over process-local source partitions, then
merges them into `symmodel`.
"""
function compute_abstract_system_distributed!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    procs = Distributed.workers(),
    nparts::Int = length(procs),
    partition_strategy::Symbol = :roundrobin,
    threaded_per_worker::Bool = false,
    verbose::Bool = true,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
)
    transitions = collect_abstract_transitions_distributed(
        symmodel,
        concrete_system_approx;
        procs = procs,
        nparts = nparts,
        partition_strategy = partition_strategy,
        threaded_per_worker = threaded_per_worker,
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
    )

    isempty(transitions) || add_transitions!(symmodel, transitions)
    return symmodel
end

