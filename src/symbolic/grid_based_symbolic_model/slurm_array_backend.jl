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

function slurm_chunk_file(outdir::AbstractString, chunk_id::Int, nchunks::Int)
    return joinpath(outdir, "transitions_part_$(chunk_id)_of_$(nchunks).bin")
end

function collect_abstract_transitions_chunk(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    chunk_id::Int,
    nchunks::Int;
    partition_strategy::Symbol = :contiguous,
    print_level::Int = 0,
    kwargs...,
)
    1 <= chunk_id <= nchunks ||
        error("Invalid chunk_id=$chunk_id. Expected value in 1:$nchunks.")

    parts = partition_source_state_ids(symmodel, nchunks; strategy = partition_strategy)
    state_ids = parts[chunk_id]

    if print_level >= 1
        ninputs = length(collect(enum_inputs(symmodel)))
        chunk_work = length(state_ids) * ninputs
        max_work = maximum(length(ids) * ninputs for ids in parts)
        min_work = minimum(length(ids) * ninputs for ids in parts)

        @info(
            "SLURM-array chunk workload",
            chunk_id = chunk_id,
            nchunks = nchunks,
            n_source_states = length(state_ids),
            ninputs = ninputs,
            state_input_checks = chunk_work,
            min_chunk_state_input_checks = min_work,
            max_chunk_state_input_checks = max_work,
        )
    end

    local_symmodel = local_symmodel_from_state_ids(symmodel, state_ids)

    transitions, metadata_pairs = collect_abstract_transitions(
        local_symmodel,
        concrete_system_approx;
        print_level = print_level,
        kwargs...,
    )

    return transitions, metadata_pairs, state_ids
end

function write_transition_chunk(
    outdir::AbstractString,
    chunk_id::Int,
    nchunks::Int,
    transitions::Vector{TransitionKey},
    metadata_pairs::Vector{Pair{TransitionKey, Any}},
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
        metadata_pairs = metadata_pairs,
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
    nchunks = execution_backend.nchunks
    chunk_id = execution_backend.chunk_id

    chunk_id === nothing && error(
        "SlurmArrayBackend.chunk_id is nothing. Provide chunk_id or read SLURM_ARRAY_TASK_ID before constructing the backend.",
    )

    if print_level >= 1
        @info(
            "Starting SLURM-array abstraction chunk",
            chunk_id = chunk_id,
            nchunks = nchunks,
            outdir = execution_backend.outdir,
        )
    end

    transitions, metadata_pairs, state_ids = collect_abstract_transitions_chunk(
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
        metadata_pairs,
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
            nmetadata = length(metadata_pairs),
            file = file,
        )
    end

    if !execution_backend.write_only
        isempty(transitions) || add_transitions!(symmodel, transitions)
        add_metadata_pairs!(symmodel, metadata_pairs)
    end

    return symmodel
end

function merge_transition_chunks!(
    symmodel::GridBasedSymbolicModel,
    outdir::AbstractString;
    nchunks::Int,
    print_level::Int = 0,
)
    all_transitions = TransitionKey[]
    all_metadata_pairs = Pair{TransitionKey, Any}[]
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

        if haskey(pairs(data), :metadata_pairs)
            append!(all_metadata_pairs, data.metadata_pairs)
        end

        total_source_states += data.n_source_states

        if print_level >= 2
            @info(
                "Read transition chunk",
                chunk_id = chunk_id,
                ntransitions = data.n_transitions,
                nmetadata =
                    haskey(pairs(data), :metadata_pairs) ? length(data.metadata_pairs) : 0,
                n_source_states = data.n_source_states,
            )
        end
    end

    isempty(all_transitions) || add_transitions!(symmodel, all_transitions)
    add_metadata_pairs!(symmodel, all_metadata_pairs)

    if print_level >= 1
        @info(
            "Merged transition chunks",
            nchunks = nchunks,
            total_source_states = total_source_states,
            total_transitions = length(all_transitions),
            total_metadata = length(all_metadata_pairs),
        )
    end

    return symmodel
end
