# ==============================================================================
# One chunk of the abstraction of the normal-form biped experiment (Slurm array task)
# ==============================================================================
#
# Runs `examples/BipedRobot/biped_4d_normal_form.jl` up to the abstraction, with the
# `SY.SlurmArrayBackend`: this task computes the transitions of its own share of the state cells
# (`SLURM_ARRAY_TASK_ID` out of `SLURM_ARRAY_TASK_COUNT`) and writes them to a file of
# `DIONYSOS_TRANSITION_OUTDIR`. `merge_chunks.jl` then gathers the files. The grids, the state box and
# every other parameter are those of the example: there is a single place to edit them.
#
# Environment:
#   DIONYSOS_TRANSITION_OUTDIR   directory of the chunk files (default: slurm_array/out/transitions)
#   DIONYSOS_NCHUNKS             number of chunks (default: SLURM_ARRAY_TASK_COUNT, the size of the array)
#   SLURM_ARRAY_TASK_ID          chunk, from 1 to the number of chunks (outside Slurm: DIONYSOS_CHUNK_ID)

import Pkg
Pkg.instantiate()

import Dionysos
const SY = Dionysos.Symbolic

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

# DIONYSOS_NCHUNKS first: re-submitting a few tasks (`--array=17`) changes SLURM_ARRAY_TASK_COUNT.
nchunks = parse(Int, get(ENV, "DIONYSOS_NCHUNKS", get(ENV, "SLURM_ARRAY_TASK_COUNT", "1")))
chunk_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", get(ENV, "DIONYSOS_CHUNK_ID", "1")))
1 <= chunk_id <= nchunks ||
    error("chunk $chunk_id of $nchunks: submit the array as --array=1-$nchunks")
outdir = abspath(
    get(ENV, "DIONYSOS_TRANSITION_OUTDIR", joinpath(@__DIR__, "out", "transitions")),
)

# Read by the example before it defaults them.
EXECUTION_BACKEND = SY.SlurmArrayBackend(
    nchunks,
    chunk_id,
    outdir,
    :contiguous,
    true, # write_only: the transitions go to the file, not to the local abstraction
)
ABSTRACTION_ONLY = true

@info "Abstraction chunk" chunk_id nchunks outdir
t = @elapsed include(EXAMPLE)
file = SY.slurm_chunk_file(outdir, chunk_id, nchunks)
@info "Chunk written" file seconds = round(t; digits = 1) size_MB =
    round(filesize(file) / 1e6; digits = 1)
