using Distributed
addprocs(4)

@everywhere using Dionysos
@everywhere using JuMP
@everywhere const DI = Dionysos
@everywhere const MP = Dionysos.Mapping
@everywhere const SY = Dionysos.Symbolic
@everywhere const ST = Dionysos.System

@everywhere include("path/to/your/distributed_abstraction.jl")


# This is the easiest partition. Later we could partition geometrically.
function partition_states(Xset_global::MP.AbstractStateSet, XMapping::MP.AbstractMapping, nparts::Int)
    states = collect(MP.enum_states(Xset_global, XMapping))
    chunks = [Int[] for _ in 1:nparts]

    for (k, q) in enumerate(states)
        push!(chunks[mod1(k, nparts)], q)
    end

    return chunks
end


@everywhere function run_distributed_job(job)
    local_sym = job.local_sym
    trans = Tuple{Int,Int,Int}[]
    compute_local_transitions!(
        trans,
        local_sym,
        job.concrete_system_approx;
        threaded = job.threaded,
        verbose = job.verbose,
        update_interval = job.update_interval,
    )
    return trans
end

function distributed_compute_abstract_system!(
    abstract_system::SY.GridBasedSymbolicModel,
    concrete_system_approx;
    workers_list = workers(),
    threaded_per_worker::Bool = false,
    verbose::Bool = true,
    update_interval::Int = Int(1e5),
)
    Xmap = SY.get_state_mapping(abstract_system)
    Xset = get_source_domain(abstract_system)
    Rset = get_retained_domain(abstract_system)

    parts = partition_states(Xset, Xmap, length(workers_list))

    jobs = map(parts) do states_part
        Xset_local = ExplicitStateSubset(states_part)

        local_sym = make_local_symbolic_model(abstract_system, Xset_local)

        (
            local_sym = local_sym,
            concrete_system_approx = concrete_system_approx,
            threaded = threaded_per_worker,
            verbose = false,
            update_interval = update_interval,
        )
    end

    results = pmap(run_distributed_job, jobs)

    for trans in results
        add_transitions!(abstract_system, trans)
    end

    return abstract_system
end

function collect_transition_tuples(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    kwargs...
)::Vector{Tuple{Int,Int,Int}}
    trans = Tuple{Int,Int,Int}[]
    _collect_transition_tuples!(trans, symmodel, concrete_system_approx; kwargs...)
    return trans
end