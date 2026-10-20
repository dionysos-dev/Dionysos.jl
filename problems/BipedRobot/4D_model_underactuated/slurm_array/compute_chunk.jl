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

# A task does not build the precompilation cache: with 300 tasks at once it would take all of them an hour,
# waiting for each other's lock files. If the cache is not valid on this node (another kind of CPU than the one
# `precompile.jl` ran on, or it did not run), stop at once, and say so. DIONYSOS_REQUIRE_CACHE=0 disables it.
# This is the first thing the task does, before `Pkg` is even loaded: in a depot without a cache (a CPU type
# nobody has precompiled on), loading `Pkg` compiles it and its standard libraries first, minutes. Nor does
# the task call `Pkg.instantiate()`: the packages are installed (by `precompile.jl`, on the shared home).
include(joinpath(@__DIR__, "cache_check.jl"))
if get(ENV, "DIONYSOS_REQUIRE_CACHE", "1") != "0"
    stale = stale_packages()
    if !isempty(stale)
        @error "The precompilation cache is not valid on this node: this task stops instead of compiling" host =
            gethostname() cpu = Sys.CPU_NAME cache_key = get(ENV, "DIONYSOS_CACHE_KEY", "-") stale hint = "run submit_precompile.sbatch on this kind of node, or keep the jobs on one node type (DIONYSOS_SBATCH_ARGS of submit_all.sh): README §6"
        exit(3)
    end
end

# The packages of the experiment are loaded (and timed) here, with `import` so that no name is brought into
# `Main`: with a valid precompilation cache this takes seconds. If it takes minutes, the cache is not valid
# on this node (see precompile.jl and README §6).
t_load = time()
import Dionysos
import StaticArrays, MathematicalSystems, LazySets, MathOptInterface, Symbolics
import IntervalArithmetic, JLD2, Plots
const SY = Dionysos.Symbolic
@info "Start-up of the chunk task" host = gethostname() cpu = Sys.CPU_NAME julia = VERSION packages_loaded_in_s =
    round(time() - t_load; digits = 1)

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
