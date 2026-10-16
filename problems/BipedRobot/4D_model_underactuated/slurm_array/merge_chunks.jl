# ==============================================================================
# Merge the chunks of the abstraction of the normal-form biped experiment
# ==============================================================================
#
# Reads the `nchunks` files that the tasks of `compute_chunk.jl` wrote in `DIONYSOS_TRANSITION_OUTDIR`
# and writes the whole abstraction to `DIONYSOS_ABSTRACTION_FILE` (`ABSTRACTION_FILE` of the example):
# run the example with that variable on any machine to load it instead of recomputing it and to carry
# on with the control synthesis and the animation. The parameters are those of the example, as for the
# chunks: they must not change in between (the file stores them and the example checks them).
#
# Environment:
#   DIONYSOS_NCHUNKS             number of chunks (the size of the array)
#   DIONYSOS_TRANSITION_OUTDIR   directory of the chunk files (default: slurm_array/out/transitions)
#   DIONYSOS_ABSTRACTION_FILE    output file (default: slurm_array/out/normal_form_abstraction.jld2)

import Pkg
Pkg.instantiate()

using Printf

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

haskey(ENV, "DIONYSOS_NCHUNKS") ||
    error("set DIONYSOS_NCHUNKS to the number of chunks of the array")
nchunks = parse(Int, ENV["DIONYSOS_NCHUNKS"])
outdir = abspath(
    get(ENV, "DIONYSOS_TRANSITION_OUTDIR", joinpath(@__DIR__, "out", "transitions")),
)

# Read by the example before it defaults them.
ABSTRACTION_FILE = abspath(
    get(
        ENV,
        "DIONYSOS_ABSTRACTION_FILE",
        joinpath(@__DIR__, "out", "normal_form_abstraction.jld2"),
    ),
)
MERGE_CHUNKS = (; outdir, nchunks)
ABSTRACTION_ONLY = true

isfile(ABSTRACTION_FILE) &&
    @warn "$ABSTRACTION_FILE exists: it is loaded and nothing is merged. Delete it to merge again."
@info "Merging the abstraction chunks" nchunks outdir ABSTRACTION_FILE
t = @elapsed include(EXAMPLE)
@printf(
    "Abstraction ready in %.1f s: %s (%.2f GB)\n",
    t,
    ABSTRACTION_FILE,
    filesize(ABSTRACTION_FILE) / 1e9
)
