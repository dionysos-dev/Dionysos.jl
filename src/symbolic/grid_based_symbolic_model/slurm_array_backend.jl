# ------------------------------------------------------------
# SLURM Array Execution
# ------------------------------------------------------------

import Serialization

"""
    SlurmArrayBackend(
        nchunks,
        chunk_id=nothing,
        outdir,
        partition_strategy=:contiguous,
        write_only=true,
    )

Execution using SLURM array jobs (file-based parallelism).

# Parameters
- `nchunks`: total number of chunks (slurm array size).
- `chunk_id`: current chunk (defaults to `SLURM_ARRAY_TASK_ID`).
- `outdir`: directory where results are written.
- `partition_strategy`: how to split states (`:contiguous` or `:roundrobin`).
- `write_only`: if `true`, only writes transitions to disk (recommended for SLURM).
"""
struct SlurmArrayBackend <: AbstractExecutionBackend
    nchunks::Int
    chunk_id::Union{Nothing, Int}
    outdir::String
    partition_strategy::Symbol
    write_only::Bool
end

SlurmArrayBackend(nchunks::Int, outdir::String) =
    SlurmArrayBackend(nchunks, nothing, outdir, :contiguous, true)

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    execution_backend::SlurmArrayBackend;
    kwargs...,
)
    return compute_abstract_system_slurm_array!(
        symmodel,
        concrete_system_approx,
        execution_backend;
        kwargs...,
    )
end

current_slurm_chunk_id() = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

current_slurm_nchunks() = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))

function slurm_chunk_file(outdir::AbstractString, chunk_id::Int, nchunks::Int)
    return joinpath(outdir, "transitions_part_$(chunk_id)_of_$(nchunks).bin")
end

function collect_abstract_transitions_chunk(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    chunk_id::Int,
    nchunks::Int;
    partition_strategy::Symbol = :contiguous,
    kwargs...,
)
    1 <= chunk_id <= nchunks ||
        error("Invalid chunk_id=$chunk_id. Expected value in 1:$nchunks.")

    parts = partition_source_state_ids(symmodel, nchunks; strategy = partition_strategy)

    state_ids = parts[chunk_id]
    local_symmodel = local_symmodel_from_state_ids(symmodel, state_ids)

    transitions =
        collect_abstract_transitions(local_symmodel, concrete_system_approx; kwargs...)

    return transitions, state_ids
end

function write_transition_chunk(
    outdir::AbstractString,
    chunk_id::Int,
    nchunks::Int,
    transitions::Vector{Tuple{Int, Int, Int}},
    state_ids::Vector{Int};
    nstates::Int,
    ninputs::Int,
)
    mkpath(outdir)

    file = slurm_chunk_file(outdir, chunk_id, nchunks)

    data = (
        chunk_id = chunk_id,
        nchunks = nchunks,
        nstates = nstates,
        ninputs = ninputs,
        n_source_states = length(state_ids),
        source_state_ids = state_ids,
        n_transitions = length(transitions),
        transitions = transitions,
    )

    open(file, "w") do io
        return Serialization.serialize(io, data)
    end

    return file
end

function read_transition_chunk(outdir::AbstractString, chunk_id::Int, nchunks::Int)
    file = slurm_chunk_file(outdir, chunk_id, nchunks)

    isfile(file) || error("Missing transition chunk file: $file")

    return open(file, "r") do io
        return Serialization.deserialize(io)
    end
end

function compute_abstract_system_slurm_array!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    execution_backend::SlurmArrayBackend;
    print_level::Int = 0,
    kwargs...,
)
    nchunks =
        execution_backend.nchunks > 0 ? execution_backend.nchunks : current_slurm_nchunks()
    chunk_id =
        execution_backend.chunk_id === nothing ? current_slurm_chunk_id() :
        execution_backend.chunk_id

    if print_level >= 1
        @info(
            "Starting SLURM-array abstraction chunk",
            chunk_id = chunk_id,
            nchunks = nchunks,
            outdir = execution_backend.outdir,
        )
    end

    transitions, state_ids = collect_abstract_transitions_chunk(
        symmodel,
        concrete_system_approx,
        chunk_id,
        nchunks;
        partition_strategy = execution_backend.partition_strategy,
        print_level = print_level,
        kwargs...,
    )

    file = write_transition_chunk(
        execution_backend.outdir,
        chunk_id,
        nchunks,
        transitions,
        state_ids;
        nstates = get_n_state(symmodel),
        ninputs = get_n_input(symmodel),
    )

    if print_level >= 1
        @info(
            "Finished SLURM-array abstraction chunk",
            chunk_id = chunk_id,
            nchunks = nchunks,
            ntransitions = length(transitions),
            file = file,
        )
    end

    if !execution_backend.write_only
        isempty(transitions) || add_transitions!(symmodel, transitions)
    end

    return symmodel
end

function merge_transition_chunks!(
    symmodel::GridBasedSymbolicModel,
    outdir::AbstractString;
    nchunks::Int,
    print_level::Int = 0,
)
    all_transitions = Tuple{Int, Int, Int}[]
    total_source_states = 0

    expected_nstates = get_n_state(symmodel)
    expected_ninputs = get_n_input(symmodel)

    for chunk_id in 1:nchunks
        data = read_transition_chunk(outdir, chunk_id, nchunks)

        data.chunk_id == chunk_id || error(
            "Chunk metadata mismatch: expected chunk_id=$chunk_id, got $(data.chunk_id).",
        )

        data.nchunks == nchunks ||
            error("Chunk $chunk_id has nchunks=$(data.nchunks), expected $nchunks.")

        data.nstates == expected_nstates || error(
            "Chunk $chunk_id has nstates=$(data.nstates), expected $expected_nstates.",
        )

        data.ninputs == expected_ninputs || error(
            "Chunk $chunk_id has ninputs=$(data.ninputs), expected $expected_ninputs.",
        )

        append!(all_transitions, data.transitions)
        total_source_states += data.n_source_states

        if print_level >= 2
            @info(
                "Read transition chunk",
                chunk_id = chunk_id,
                ntransitions = data.n_transitions,
                n_source_states = data.n_source_states,
            )
        end
    end

    isempty(all_transitions) || add_transitions!(symmodel, all_transitions)

    if print_level >= 1
        @info(
            "Merged transition chunks",
            nchunks = nchunks,
            total_source_states = total_source_states,
            total_transitions = length(all_transitions),
        )
    end

    return symmodel
end
