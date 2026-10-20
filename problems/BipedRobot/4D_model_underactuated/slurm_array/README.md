# Slurm-array abstraction

The abstraction of the normal-form biped experiment
([`examples/BipedRobot/biped_4d_normal_form.jl`](../../../../examples/BipedRobot/biped_4d_normal_form.jl))
is by far its most expensive part: one image to compute per `(cell, input)` pair (27 inputs by default, `(2k+1)³` with `input_steps = k`), and a
number of cells that is the product of the cell counts along the five axes. Here it is split over the
tasks of a **Slurm job array** (`SY.SlurmArrayBackend`), one core each, as in
`problems/BipedRobot/execution/slurm_array`; the transitions are then merged into one file, and the
control synthesis and the animation then run on the server too (so that only the animation has to be
downloaded), from the saved file, which can be reused to rerun the synthesis with other starts or tolerances.

```
submit_all.sh N ──► precompile.jl ──► array 1…N : compute_chunk.jl ──► N chunk files ──► merge_chunks.jl ──► .jld2 ──► solve.jl ──► .gif
                    (compute node,    (1 core each)                     (outdir)          (saved, reusable)
                     once)             each step runs only if the previous one succeeded
```

| file | role |
| :--- | :--- |
| [`submit_all.sh`](submit_all.sh) | submits the precompile job, the array of `N` chunks, the merge and the solve job, each after the previous one has succeeded (`N` is chosen here only) |
| [`submit_precompile.sbatch`](submit_precompile.sbatch) / [`precompile.jl`](precompile.jl) | builds the Julia precompilation cache on a compute node, once, and checks that it is valid there |
| [`env.sh`](env.sh), [`cache_check.jl`](cache_check.jl) | where the cache lives (one directory per CPU type) and the check that the tasks make before starting |
| [`status.sh`](status.sh) | what happened to the jobs of the last `submit_all.sh`: state, exit code, peak memory, end of the logs (§6) |
| [`submit_compute_chunks.sbatch`](submit_compute_chunks.sbatch) / [`compute_chunk.jl`](compute_chunk.jl) | one task of the array: computes its share of the source cells and writes `transitions_part_<id>_of_<N>.bin` |
| [`submit_merge_chunks.sbatch`](submit_merge_chunks.sbatch) / [`merge_chunks.jl`](merge_chunks.jl) | copies the `N` chunks, one after the other, into the abstraction file `DIONYSOS_ABSTRACTION_FILE` (a few GB of memory whatever the size) |
| [`submit_solve.sbatch`](submit_solve.sbatch) / [`solve.jl`](solve.jl) | control synthesis, closed loop, figures and animation on the saved abstraction; rerunnable alone |
| `logs/` | job logs (not versioned) |
| `out/` | default location of the chunks, of the merged file and of the results in `out/results/` (`*.jld2` and `out/` are ignored by git) |

There is no copy of the parameters: **the grids, the state box, `tstep`, the number of inputs and the
obstacle are those of the example**, so edit the example (or set the variables the example reads, see
below) and nothing else. They are stored in the merged file, and the example refuses to load a file
computed for other parameters. **Do not change them between the chunks and the merge** (nor between two
tasks of the array): the chunks only record the number of states and of inputs.

## 1. Prepare the environment (once, on the login node)

```
julia --project=problems/BipedRobot -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

This installs the packages (compute nodes often have no internet access) and downloads what the
animation needs (Plots, the GIF encoder): the solve job writes it headless (`GKSwstype=100`).

**The precompilation cache is not portable between machines**: it holds machine code and is rejected
on another CPU model (or Julia, or `JULIA_CPU_TARGET`), then rebuilt. Precompiling on the login node
therefore says nothing about the compute nodes. `submit_all.sh` submits `submit_precompile.sbatch` first
and the array waits for it: it builds the cache on a compute node, for the packages of the experiment
only (not the whole project), then loads them in a fresh process, twice, and fails (`PRECOMPILE FAILED`)
if the second load is slow, which means the cache is rejected on that node. Without it, each of the
`N` tasks would compile for itself, minutes of its budget, all of them fighting for the same lock files (the
`attempting to remove probably stale pidfile` warnings). The tasks themselves never precompile the whole
project: the merge and solve jobs call `Pkg.instantiate()` with `JULIA_PKG_PRECOMPILE_AUTO=0`, and the tasks of the array do not call it (nor load `Pkg`). `DIONYSOS_SKIP_PRECOMPILE=1`
skips the job (when the cache is known to be valid).

**Heterogeneous clusters.** A cache file has the same name whatever the CPU that built it. On a cluster
whose nodes are not all alike (check: `sinfo -p <partition> -o "%20N %f"`, or the node names in `squeue`) and
that shares the home directory, jobs on different CPU types overwrite each other's cache: each finds the
cache of the other invalid and recompiles, and they all wait for each other's lock files (`Being precompiled
by another machine`). That is what makes tasks hang for an hour (§6). The pipeline protects itself:

* the jobs use a depot of their own for the cache, **one directory per CPU type** (`env.sh`, in
  `~/.julia_biped_cache/<cpu>/`, in front of `~/.julia` where the sources of the packages stay): jobs on the
  same type share it, those on other types have theirs, nothing is overwritten;
* every task of the array **checks the cache before doing anything** and, if it is not valid on its node,
  stops at once with exit code 3 and says so, instead of compiling for an hour (`DIONYSOS_REQUIRE_CACHE=0`
  disables it);
* the precompile job builds the cache for the CPU type of the node it runs on, so **it and the array must be
  on the same type**: keep the whole pipeline on one, with the extra `sbatch` arguments of `submit_all.sh`
  (applied to every job):

```
DIONYSOS_SBATCH_ARGS="--partition=cp3 --constraint=Skylake" \
    problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 300
```

  (`--constraint` takes a feature of `sinfo -o "%N %f"`; without features use `--exclude=<nodes>` or
  `--nodelist=<nodes>`; on the UCLouvain cluster the features of `cp3` are `SkyLake`, `CascadeLake`,
  `Rome` and `Genoa`). `submit_all.sh` **refuses to submit** when the cluster has nodes with different
  features and nothing pins the jobs (nor a portable `JULIA_CPU_TARGET`); `DIONYSOS_ANY_NODE=1` overrides it;
* or, to use all the node types, export one portable CPU target before `submit_all.sh` (sbatch passes it
  on); the cache is then one for all the types (directory `portable`), and the tasks that find it invalid on
  their node stop at once as above. I have not tested this target on your hardware:

```
export JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"   # the target of the official Julia builds
```

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
threads). A task needs, besides its share of the work, a fixed **start-up** (loading Julia and the packages, building
the model and the local growth bound, carving the state box: of the order of a minute with a valid
precompilation cache, §1; the chunk log says `packages_loaded_in_s`), so the
chunks must be long enough for it not to dominate: aim for 10–30 minutes per task, i.e.
`N ≈ total core time / 20 min`. Set `--time` of `submit_compute_chunks.sbatch` (default 1 h, 4 GB) with a
margin of 2×; the chunks have the same number of source cells, but the work per cell varies a little.

**Memory.** The three jobs do not need the same:

* the **chunk** tasks need little (4 GB by default): their share of the transitions, 1/N of the total;
* the **merge** copies the chunks one after the other into the file, without building the abstraction:
  the memory is that of the model (about 2 GB) plus one chunk (about 24 bytes per transition *of that chunk*),
  so a few GB whatever the size of the abstraction. (It used to build the whole automaton, hundreds of bytes
  per transition: with billions of transitions that is hundreds of GB and the job was killed. Estimate the
  transitions from the chunks, `du -sh out/transitions`, 23 bytes each);
* the **solve** job has to put the whole abstraction into an automaton, as any synthesis does: about
  **150 bytes per transition** (measured: a peak of 2.8 GB for 9.3 M transitions, of which about 1.2 GB is
  the model, synthesis and animation included). With `T` transitions in the file (the merge log says it:
  `MERGE OK … T transitions`), ask for about `2 GB + 200 bytes × T`, i.e. 400 GB for 2 billion
  transitions. This is the one that decides whether the abstraction fits the cluster nodes.

The merged file takes 12 bytes per transition. A grid twice as fine along the ankle angle and momentum
has 4 times more cells and transitions.

## 3. Submit

From the root of the repository, in a shell where `julia` works (the jobs inherit its `PATH`; over `ssh host 'command'` use a login shell, `ssh host 'bash -l -c "..."'`; the scripts stop at once with a message if `julia` is missing):

```
problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200        # N = 200 chunks
problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200 50     # at most 50 tasks at once
```

This submits four jobs, each one running only if the previous one has succeeded (`--dependency=afterok`):
the precompilation (§1), the array `1-200`, the merge (`DIONYSOS_NCHUNKS=200`) and the **solve** job (control synthesis, closed
loop, figures and animation, see §4), and prints the name of the animation it will write
(`slurm_array/out/results/biped_4d_normal_form_job<solve jobid>.gif`: that is the file to download).
`DIONYSOS_SKIP_SOLVE=1` submits up to the merge only. If a job fails, the ones after it are
cancelled (`--kill-on-invalid-dep`): `status.sh` (§6) tells which one failed and why.

**Only the merge** (the chunks are all there): yes, `sbatch` it, **from the root of the repository** (the log
paths of the script are relative to it, and the `logs/` folder must exist):

```
sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_merge_chunks.sbatch
```

The number of chunks is read from the names of the chunk files (`DIONYSOS_NCHUNKS=300` overrides it), and
the merge stops at once, naming them, if some are missing. Set `DIONYSOS_TRANSITION_OUTDIR` and
`DIONYSOS_ABSTRACTION_FILE` again if you had changed them for `submit_all.sh`. If the output file exists, the
job loads it (loading is the "merge" then, which is not what you want): delete it, and the `….jld2.tmp` of a
merge that was killed. The same goes for **only the solve job** (§4), and to re-run some failed chunk tasks
only: `DIONYSOS_NCHUNKS=300 sbatch --array=17,42 …/submit_compute_chunks.sbatch` (the chunk files of the other
tasks are kept), then the merge.

Directories: the chunks go to `slurm_array/out/transitions/`, the merged file to
`slurm_array/out/normal_form_abstraction.jld2` and the results to `slurm_array/out/results/`. Choose others (e.g. a scratch directory) with the
environment of the submission:

```
DIONYSOS_TRANSITION_OUTDIR=$GLOBALSCRATCH/biped4d/transitions \
DIONYSOS_ABSTRACTION_FILE=$GLOBALSCRATCH/biped4d/abstraction.jld2 \
DIONYSOS_OUTPUT_DIR=$HOME/biped4d_results \
    problems/BipedRobot/4D_model_underactuated/slurm_array/submit_all.sh 200
```

`submit_all.sh` refuses to start if the chunk directory already has chunk files, which would be
merged silently with the new ones; delete them (they are only useful for the run that made them). If
the merged file exists, the merge job loads it instead of merging: delete it to merge again.

## 4. The solve job, and rerunning it

`submit_solve.sbatch` runs the example on the merged file (`solve.jl`): it loads the transitions (nothing
is computed; if the file is missing the job stops instead of computing the abstraction sequentially),
synthesizes the controller, simulates the closed loop and writes, in `DIONYSOS_OUTPUT_DIR`,
`biped_4d_normal_form_<scenario>.gif` (the animation), `….png` (snapshots) and `…_ankle.png`, with the
log in `logs/solve_<jobid>.out` (`footstep success`, number of steps, ankle range). If there is no
controller from the start state, the job says so and ends with exit code 2 (no gif).

**Memory.** The solve job puts the whole abstraction in an automaton: see the estimate of §2
(`--mem` of `submit_solve.sbatch`, 250 GB as a placeholder; 2 CPUs and 4 h too). Its log says how many
transitions the file has and loads them part by part with the peak memory so far, so a job that is
killed leaves the point it had reached.

**Rerun with another start or tolerance.** What does not enter the abstraction can change freely: the
foothold tolerance, the initial configuration and momentum, the target momentum (and the foothold). Put
the new values in a Julia file,

```julia
# my_params.jl
FOOTHOLD_TOLERANCE = 0.03
INITIAL_VELOCITY = (-0.44, -0.36)
Q0 = SVector(0.0, 0.5, -0.4, -0.9)
TARGET_VELOCITY = (-0.3, -0.1)
```

and submit the solve job alone, from the root of the repository:

```
DIONYSOS_SOLVE_PARAMS=my_params.jl DIONYSOS_SCENARIO=tol3cm \
    sbatch problems/BipedRobot/4D_model_underactuated/slurm_array/submit_solve.sbatch
```

`DIONYSOS_SCENARIO` names the outputs (default `job<jobid>`), so that the reruns do not overwrite each
other; the params file is copied next to them as `…_params.jl`. Submit several at once to scan values.
`DIONYSOS_ABSTRACTION_FILE` and `DIONYSOS_OUTPUT_DIR` are the variables of the merge, if you changed them.

What **does** enter the abstraction (and is checked against the file: it is refused if it differs) is the
example's `DISC`, `X_BOUNDS`, `ONE_JOINT_PER_STEP` and obstacle: change them in the example and a new
abstraction is needed. A different `X_BOUNDS` on the same grid is one too: compute one big box once and
try several starts in it.

The file holds the transitions only (`target`, `source`, `symbol` as integer vectors, one part per chunk,
plus the parameters): the compiled symbolic dynamics cannot be saved, they are rebuilt from the model at each
run (a few seconds). Files written by the first version (one part) are still read.

### Locally

On a machine that can hold the abstraction, copy the merged file (`scp`/`rsync`) and run the example with
the same parameters, pointing it at it (the file must be an absolute path):

```julia
ABSTRACTION_FILE = "/path/to/normal_form_abstraction.jld2"
include("examples/BipedRobot/biped_4d_normal_form.jl")
```

The example prints `Loading the abstraction from … (not computed)`.

## 5. How the pieces fit

The variables of the example that these scripts set before including it:

| variable | chunk task | merge |
| :--- | :--- | :--- |
| `EXECUTION_BACKEND` | `SY.SlurmArrayBackend(N, chunk_id, outdir, :contiguous, true)` (the task writes its chunk and keeps nothing) | default |
| `MERGE_CHUNKS` | – | `(; outdir, N)`: the transitions are read from the chunks instead of computed, and written to `ABSTRACTION_FILE` (the run has to stop there: `ABSTRACTION_ONLY`) |
| `ABSTRACTION_FILE` | – | where the merged abstraction is written (the solve job: where it is read) |
| `ABSTRACTION_ONLY` | `true` | `true` |
| `OUTPUT_DIR`, `SCENARIO_NAME`, `RUN_PLOTS` | – | – (set by `solve.jl`: where and under which name the figures are written) |

They are also what you would set to run the steps by hand, e.g. a chunk outside Slurm
(`DIONYSOS_CHUNK_ID=2 DIONYSOS_NCHUNKS=3 julia --project=problems/BipedRobot .../compute_chunk.jl`). The
functions behind them are in `LocalGrowthAbstraction.jl`: `build_abstraction(...; merge_chunks,
abstraction_file, abstraction_info)`, `merge_chunks!` (into the optimizer), `merge_chunks_to_file` (chunks to
file, what the merge job uses), `save_abstraction` and `load_abstraction!`.

## 6. Did it work? Reading the jobs

```
problems/BipedRobot/4D_model_underactuated/slurm_array/status.sh                  # the last submit_all.sh
problems/BipedRobot/4D_model_underactuated/slurm_array/status.sh CHUNKS MERGE [SOLVE]   # given job ids
```

shows the chunk tasks by state (and the ones that did not complete), the number of chunk files, and for
the merge and solve jobs the `sacct` line (state, exit code, elapsed, requested and **peak memory**) and the
end of their logs. The same by hand: `sacct -j <jobid> --format=JobID,JobName,State,ExitCode,Elapsed,ReqMem,MaxRSS`.

| `State` | meaning |
| :--- | :--- |
| `COMPLETED` | the job ended with status 0 |
| `OUT_OF_MEMORY` | it needed more than `--mem`: Slurm killed it, and a killed job logs nothing, its log just stops (the last `Merging chunk k / N … peak_memory_GB` or `Loaded part k` line says how far it got). `MaxRSS` is the peak |
| `TIMEOUT` | it needed more than `--time` |
| `FAILED` | the script ended with an error: the `.err` log has it |
| `CANCELLED` | by you, or because the job it waited for failed: look at the first job that is not `COMPLETED` |

The merge and solve scripts log each stage with a time stamp and end with one of
`MERGE OK … transitions, … peak memory`, `MERGE FAILED …`, `SOLVE OK: <gif>`, `SOLVE FAILED …`; the sbatch
scripts add `[merge] julia exited with status …`. In the **merge**, the usual failures are chunk files that are
missing (named at once, before any work), chunks computed for another grid or state box (the number of states
differs), and a merged file that already exists (it is loaded, not merged: a warning says so). In the
**solve** job: the abstraction file refused because the example's parameters are not the ones it was
computed for (the message prints both), no controller from the start state (exit code 2: raise the
tolerance, change the start), or `OUT_OF_MEMORY` while loading.

### The tasks take much longer than their share of the work, or hang

Look at a chunk log (`logs/chunk_<jobid>_<task>.err`/`.out`) and at `squeue` (which nodes?):

* `Precompiling packages…`, `Being precompiled by another process / machine` for packages of the
  experiment (even JuMP, Polyhedra, …) and nothing else for tens of minutes: the tasks are compiling the
  packages themselves, all at once, waiting for each other's lock files. If the tasks run on **several node
  families** (`mb-sky…`, `mb-cas…`, `mb-rom…`) the cache is being rebuilt back and forth by each CPU type:
  cancel the jobs (`scancel <array job> <merge job> <solve job>`: they would only reach the time limit) and
  keep the pipeline on one node type or use a portable target (§1). The current scripts do not do this any
  more: the cache is per CPU type, and a task without a valid cache stops (exit 3, `The precompilation cache
  is not valid on this node`).
* `The precompilation cache is not valid on this node`: the precompile job ran on another CPU type (compare
  `CPU …` in `logs/precompile_<jobid>.out` with `cpu =` in the chunk log), or it did not run or fail
  (`PRECOMPILE OK` or `PRECOMPILE FAILED` at the end of its log). Re-run the failed tasks after building the
  cache for their type: `sbatch --constraint=<that type> …/submit_precompile.sbatch`, then
  `DIONYSOS_NCHUNKS=N sbatch --array=<failed tasks> …/submit_compute_chunks.sbatch`.
* `Start-up of the chunk task … packages_loaded_in_s = …` is the time to load the packages: seconds to a
  minute with a valid cache. The other part of the start-up is the model, compiled in every task (the `@time`
  lines of the log: a few tens of seconds locally, more on a slower node): that one is not cached, it is the
  price of each task.
