# ==============================================================================
# Merge the chunks of the abstraction of the normal-form biped experiment
# ==============================================================================
#
# Reads the chunk files that the tasks of `compute_chunk.jl` wrote in `DIONYSOS_TRANSITION_OUTDIR` and
# writes the whole abstraction to `DIONYSOS_ABSTRACTION_FILE` (`ABSTRACTION_FILE` of the example). The
# chunks are copied into the file one after the other: the abstraction is not built in memory, so the
# job needs the memory of one chunk (a few GB), not of the abstraction. Run the example with that
# file on a machine that can hold the abstraction (or `solve.jl`) to load it and carry on with the
# control synthesis and the animation. The parameters are those of the example, as for the chunks:
# they must not change in between (the file stores them and the example checks them).
#
# Everything is logged with a time stamp and the peak memory of the process; a failure ends with a
# `MERGE FAILED` line and exit code 1. (A job killed for lack of memory cannot log anything: its log
# stops abruptly, and `sacct` says `OUT_OF_MEMORY`, see README §5.)
#
# Environment:
#   DIONYSOS_TRANSITION_OUTDIR   directory of the chunk files (default: slurm_array/out/transitions)
#   DIONYSOS_NCHUNKS             number of chunks (default: read from the names of the chunk files)
#   DIONYSOS_ABSTRACTION_FILE    output file (default: slurm_array/out/normal_form_abstraction.jld2)

# `instantiate` must not precompile the whole project (MeshCat, RigidBodyDynamics, ... that the experiment
# never loads): the precompilation is the job of `precompile.jl`, once, before the array.
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
import Pkg
Pkg.instantiate()

using Printf
using Dates
import JLD2

stamp(msg) = @info "[merge $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
peak_gb() = round(Sys.maxrss() / 2^30; digits = 1)

const EXAMPLE = abspath(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "..",
        "examples",
        "BipedRobot",
        "biped_4d_normal_form.jl",
    ),
)

# The number of chunks, from the names `transitions_part_<id>_of_<N>.bin`; the missing ones are
# reported now, not after a minute of setup.
function chunk_files(outdir)
    isdir(outdir) || error("$outdir does not exist: no chunk was written there (DIONYSOS_TRANSITION_OUTDIR)")
    found = Dict{Int, Vector{Int}}()
    for name in readdir(outdir)
        m = match(r"^transitions_part_(\d+)_of_(\d+)\.bin$", name)
        m === nothing || push!(get!(found, parse(Int, m[2]), Int[]), parse(Int, m[1]))
    end
    isempty(found) && error("no chunk file transitions_part_<id>_of_<N>.bin in $outdir")
    return found
end

try
    stamp("start: job $(get(ENV, "SLURM_JOB_ID", "-")), memory requested " *
          "$(get(ENV, "SLURM_MEM_PER_NODE", "?")) MB, Julia $(VERSION), $(gethostname())")
    outdir = abspath(
        get(ENV, "DIONYSOS_TRANSITION_OUTDIR", joinpath(@__DIR__, "out", "transitions")),
    )
    found = chunk_files(outdir)
    nchunks = if haskey(ENV, "DIONYSOS_NCHUNKS")
        parse(Int, ENV["DIONYSOS_NCHUNKS"])
    elseif length(found) == 1
        only(keys(found))
    else
        error(
            "chunk files for several numbers of chunks in $outdir ($(sort(collect(keys(found))))): " *
            "remove the old ones, or set DIONYSOS_NCHUNKS",
        )
    end
    ids = sort(get(found, nchunks, Int[]))
    absent = setdiff(1:nchunks, ids)
    isempty(absent) || error(
        "$(length(absent)) of the $nchunks chunks are missing in $outdir (first ones: " *
        "$(first(absent, 20))): did all the tasks of the array succeed?",
    )
    chunk_bytes = sum(filesize(joinpath(outdir, "transitions_part_$(i)_of_$(nchunks).bin")) for i in ids)
    stamp("$nchunks chunks found in $outdir, $(round(chunk_bytes / 2^30; digits = 1)) GB " *
          "(about $(round(Int, chunk_bytes / 23)) transitions)")

    # Read by the example before it defaults them.
    global ABSTRACTION_FILE = abspath(
        get(
            ENV,
            "DIONYSOS_ABSTRACTION_FILE",
            joinpath(@__DIR__, "out", "normal_form_abstraction.jld2"),
        ),
    )
    global MERGE_CHUNKS = (; outdir, nchunks)
    global ABSTRACTION_ONLY = true

    isfile(ABSTRACTION_FILE) &&
        @warn "$ABSTRACTION_FILE exists: it is loaded and nothing is merged. Delete it to merge again."
    isfile(ABSTRACTION_FILE * ".tmp") &&
        @warn "$(ABSTRACTION_FILE).tmp is left by a merge that did not finish: it is overwritten"
    stamp("setting up the model (the merge starts after it); output $ABSTRACTION_FILE")
    t = @elapsed include(EXAMPLE)

    n = JLD2.jldopen(f -> f["n_transitions"], ABSTRACTION_FILE, "r")
    stamp(@sprintf("MERGE OK in %.1f s: %d transitions, %s (%.2f GB), peak memory %.1f GB",
        t, n, ABSTRACTION_FILE, filesize(ABSTRACTION_FILE) / 2^30, Sys.maxrss() / 2^30))
catch e
    stamp("MERGE FAILED after peak memory $(peak_gb()) GB: " * (e isa OutOfMemoryError ?
        "OutOfMemoryError: raise --mem of submit_merge_chunks.sbatch" : sprint(showerror, e)))
    @error "MERGE FAILED" exception = (e, catch_backtrace())
    exit(1)
end
