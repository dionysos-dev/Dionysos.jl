# Slurm-array abstraction

The abstraction of the normal-form biped experiment
([`examples/BipedRobot/biped_4d_normal_form.jl`](../../../../examples/BipedRobot/biped_4d_normal_form.jl))
is by far its most expensive part: one image to compute per `(cell, input)` pair (27 inputs by default, `(2k+1)³` with `input_steps = k`), and a
number of cells that is the product of the cell counts along the five axes. Here it is split over the
tasks of a **Slurm job array** (`SY.SlurmArrayBackend`), one core each, as in
`problems/BipedRobot/execution/slurm_array`; the transitions are then merged into one file, and the
control synthesis and the animation run on any machine from that file.

```
submit_all.sh N ──► array 1…N : compute_chunk.jl ──► N chunk files ──► merge_chunks.jl ──► one .jld2 ──► the example
                    (1 core each)                     (outdir)          (afterok, 1 core)                (anywhere)
```

| file | role |
| :--- | :--- |
| [`submit_all.sh`](submit_all.sh) | submits the array of `N` chunks, then the merge once all have succeeded (`N` is chosen here only) |
| [`submit_compute_chunks.sbatch`](submit_compute_chunks.sbatch) / [`compute_chunk.jl`](compute_chunk.jl) | one task of the array: computes its share of the source cells and writes `transitions_part_<id>_of_<N>.bin` |
| [`submit_merge_chunks.sbatch`](submit_merge_chunks.sbatch) / [`merge_chunks.jl`](merge_chunks.jl) | reads the `N` chunks and writes the abstraction to `DIONYSOS_ABSTRACTION_FILE` |
| `logs/` | job logs (not versioned) |
| `out/` | default location of the chunks and of the merged file (`*.jld2` and `out/` are ignored by git) |

There is no copy of the parameters: **the grids, the state box, `tstep`, the number of inputs and the
obstacle are those of the example**, so edit the example (or set the variables the example reads, see
below) and nothing else. They are stored in the merged file, and the example refuses to load a file
computed for other parameters. **Do not change them between the chunks and the merge** (nor between two
tasks of the array): the chunks only record the number of states and of inputs.

## 1. Prepare the environment (once, on the login node)

```
julia --project=problems/BipedRobot -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

The tasks run `Pkg.instantiate()` as well, but compute nodes often have no internet access, and 100
tasks precompiling at once on a shared depot would be a mess.

## 2. Choose the parameters and the number of chunks

Set the discretization and the box in the example, e.g. a finer grid:

```julia
disc = RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 1 / 60)
```

Running the example locally prints the size of the abstraction before computing it, so the run can be
interrupted right after that line:

```
julia --project=problems/BipedRobot examples/BipedRobot/biped_4d_normal_form.jl
...
State space: 123420 cells (...), 3332340 (cell, input) pairs: about 13.3 M transitions and 1.1 min on 10 thread(s)
```

Run with one thread (`-t 1`), the estimate is the total core time: 11.1 min for the default box, and the
three chunks of a local test took 3 min each (the estimate is rough, and was optimistic by 1.5× on 10
threads). A task needs, besides its share of the work, a fixed **start-up** (loading Julia and Symbolics,
building the model and the local growth bound, carving the state box: of the order of a minute), so the
chunks must be long enough for it not to dominate: aim for 10–30 minutes per task, i.e.
`N ≈ total core time / 20 min`. Set `--time` of `submit_compute_chunks.sbatch` (default 1 h, 4 GB) with a
margin of 2×; the chunks have the same number of source cells, but the work per cell varies a little.

The merge keeps all the transitions in memory: about `2 GB + 300 bytes × transitions` (measured: 3.5 GB
for 9.3 M transitions); set `--mem` of `submit_merge_chunks.sbatch` (default 16 GB) accordingly. The merged
file takes 12 bytes per transition. A grid twice as fine along the ankle angle and momentum has 4 times more
cells and transitions.

## 3. Submit

From the root of the repository:

```
problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200        # N = 200 chunks
problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200 50     # at most 50 tasks at once
```

This submits the array `1-200` and, with `--dependency=afterok`, the merge (`DIONYSOS_NCHUNKS=200`).
If a task fails, the merge does not start (`afterok`): look at `logs/chunk_<jobid>_<task>.err`, re-run
only the failed tasks with the same `N` (the chunk files of the others are kept), then submit the merge by hand:

```
DIONYSOS_NCHUNKS=200 sbatch --array=17,42 problems/BipedRobot/4D_model_underactuated/slurm_array/submit_compute_chunks.sbatch
DIONYSOS_NCHUNKS=200 sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_merge_chunks.sbatch
```

(add `--dependency=afterok:<jobid>` to the second one to wait for the first.)

Directories: the chunks go to `slurm_array/out/transitions/` and the merged file to
`slurm_array/out/normal_form_abstraction.jld2`. Choose others (e.g. a scratch directory) with the
environment of the submission:

```
DIONYSOS_TRANSITION_OUTDIR=$GLOBALSCRATCH/biped4d/transitions \
DIONYSOS_ABSTRACTION_FILE=$GLOBALSCRATCH/biped4d/abstraction.jld2 \
    problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200
```

`submit_all.sh` refuses to start if the chunk directory already has chunk files, which would be
merged silently with the new ones; delete them (they are only useful for the run that made them). If
the merged file exists, the merge job loads it instead of merging: delete it to merge again.

## 4. Use it locally

Copy the merged file (`scp`/`rsync`) and run the example with the same parameters, pointing it at it:

```julia
ABSTRACTION_FILE = "/path/to/normal_form_abstraction.jld2"
include("examples/BipedRobot/biped_4d_normal_form.jl")
```

The example then skips the computation of the abstraction (a few seconds to load) and carries on with the
start state, the control synthesis, the closed loop and the animation, as usual. **The target (`foothold`,
`FOOTHOLD_TOLERANCE`), the start (`Q0`, `INITIAL_VELOCITY`) and the target velocity can be changed
freely without recomputing the abstraction**: they do not enter the file. What does enter it (and is
checked) is `DISC`, `X_BOUNDS`, `ONE_JOINT_PER_STEP` and the obstacle. A different `X_BOUNDS` on the same grid is a different abstraction too: compute
one big box once and try several starts in it.

The file holds the transitions only (`target`, `source`, `symbol` as integer vectors, plus the
parameters): the compiled symbolic dynamics cannot be saved, they are rebuilt from the model at each run
(a few seconds).

## 5. How the pieces fit

The variables of the example that these scripts set before including it:

| variable | chunk task | merge |
| :--- | :--- | :--- |
| `EXECUTION_BACKEND` | `SY.SlurmArrayBackend(N, chunk_id, outdir, :contiguous, true)` (the task writes its chunk and keeps nothing) | default |
| `MERGE_CHUNKS` | – | `(; outdir, N)`: the transitions are read from the chunks instead of computed |
| `ABSTRACTION_FILE` | – | where the merged abstraction is written |
| `ABSTRACTION_ONLY` | `true` | `true` |

They are also what you would set to run the steps by hand, e.g. a chunk outside Slurm
(`DIONYSOS_CHUNK_ID=2 DIONYSOS_NCHUNKS=3 julia --project=problems/BipedRobot .../compute_chunk.jl`). The
functions behind them are in `LocalGrowthAbstraction.jl`: `build_abstraction(...; merge_chunks,
abstraction_file, abstraction_info)`, `merge_chunks!`, `save_abstraction` and `load_abstraction!`.
