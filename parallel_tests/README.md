# Dionysos Parallel Tests

Benchmark any Dionysos example across four parallelism strategies, collect JSON metrics per mode, and optionally generate a LaTeX/PDF comparison report.

## Parallelism Modes

| Mode | `DIONYSOS_DISTRIBUTED` | `DIONYSOS_THREADED` | What happens |
|------|----------------------|-------------------|-------------|
| **Serial** | `false` | `false` | Single thread, single process — baseline |
| **Threaded** | `false` | `true` | `Threads.@threads` on one process |
| **Distributed** | `true` | `false` | `Distributed.pmap` across workers, each sequential |
| **Hybrid** | `true` | `true` | `pmap` across workers, each using `Threads.@threads` |

## Quick Start

### Local (no SLURM)

```bash
# Run all 4 modes for the DC-DC example
bash parallel_tests/run_local.sh utils/example_distributed.jl

# Run all 4 modes with 4 distributed partitions + generate report
bash parallel_tests/run_local.sh utils/example_distributed.jl 4 --report

# Robot example
bash parallel_tests/run_local.sh problems/BipedRobot/6D_model/robot_example.jl 4 --report
```

### Supercomputer (SLURM)

```bash
# Submit all 4 modes + generate report after completion
bash parallel_tests/slurm/run_all.sh utils/example_distributed.jl --report

# Submit all 4 modes without report
bash parallel_tests/slurm/run_all.sh utils/example_distributed.jl

# Submit individual modes
sbatch parallel_tests/slurm/run_serial.sh utils/example_distributed.jl
sbatch parallel_tests/slurm/run_threaded.sh utils/example_distributed.jl
sbatch parallel_tests/slurm/run_distributed.sh utils/example_distributed.jl
sbatch parallel_tests/slurm/run_hybrid.sh utils/example_distributed.jl

# Generic: pick mode via arguments
sbatch --ntasks=1 --cpus-per-task=8 \
    parallel_tests/slurm/run_example.sh utils/example_distributed.jl false true
```

### Generate Report Separately

If you already have results and just want a report:

```bash
julia --project=. parallel_tests/src/generate_report.jl \
    parallel_tests/results/example_distributed \
    parallel_tests/report/example_distributed
```

## Command Reference

### `run_local.sh`

```
bash parallel_tests/run_local.sh <example.jl> [NPARTS] [--report]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `<example.jl>` | Yes | — | Path to Julia script (relative to project root) |
| `[NPARTS]` | No | `4` | Number of distributed partitions |
| `--report` | No | off | Generate LaTeX/PDF report after all 4 runs |

### `slurm/run_all.sh`

```
bash parallel_tests/slurm/run_all.sh <example.jl> [--report]
```

Submits 4 SLURM jobs (one per mode). With `--report`, submits a 5th dependent job that generates the report after all benchmarks complete.

### `slurm/run_example.sh` (generic)

```
sbatch [SLURM_OPTS] parallel_tests/slurm/run_example.sh <example.jl> <distributed> <threaded>
```

| Argument | Values | Description |
|----------|--------|-------------|
| `<example.jl>` | path | Julia script relative to project root |
| `<distributed>` | `true`/`false` | Enable distributed workers |
| `<threaded>` | `true`/`false` | Enable multi-threading |

Set `--ntasks` and `--cpus-per-task` via `sbatch` options to control worker/thread count.

### `src/generate_report.jl`

```
julia --project=. parallel_tests/src/generate_report.jl [results_dir] [report_dir]
```

Reads `metrics_*.json` files from `results_dir`, generates `benchmark_report.tex` (and PDF if `pdflatex` is available) in `report_dir`.

## Environment Variables

The example scripts and wrapper read these environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `DIONYSOS_DISTRIBUTED` | `false` | Enable distributed workers (`addprocs`) |
| `DIONYSOS_THREADED` | `false` | Enable multi-threading (`Threads.@threads`) |
| `DIONYSOS_NPARTS` | `8` | Number of distributed partitions / workers |
| `DIONYSOS_OUTDIR` | script dir | Output directory for results, figures, metrics JSON |

## Project Structure

```
parallel_tests/
├── src/
│   ├── run_and_collect.jl     # Wrapper: runs example + writes JSON metrics
│   └── generate_report.jl     # Reads metrics JSONs → LaTeX/PDF report
├── slurm/
│   ├── run_serial.sh          # sbatch: serial mode
│   ├── run_threaded.sh        # sbatch: threaded mode
│   ├── run_distributed.sh     # sbatch: distributed mode
│   ├── run_hybrid.sh          # sbatch: hybrid mode
│   ├── run_example.sh         # sbatch: generic (takes distributed/threaded args)
│   └── run_all.sh             # Submit all 4 modes + optional report
├── run_local.sh               # Run all 4 modes locally (no SLURM)
├── .gitignore
└── README.md
```

### Output Layout

```
parallel_tests/
├── results/<example_name>/
│   ├── serial/
│   │   ├── run.log                        # Full stdout/stderr
│   │   ├── metrics_serial_<id>.json       # Structured metrics
│   │   └── <example outputs>              # Figures, JLD2, etc.
│   ├── threaded/
│   ├── distributed/
│   └── hybrid/
├── report/<example_name>/
│   ├── benchmark_report.tex               # LaTeX source
│   └── benchmark_report.pdf               # Compiled report (if pdflatex available)
└── logs/                                  # SLURM .out/.err files
```

## Metrics JSON Format

Each run produces a `metrics_<mode>_<id>.json` file containing:

```json
{
    "mode": "hybrid",
    "example": "example_distributed.jl",
    "elapsed_sec": 5.0617,
    "wall_clock_sec": 42.31,
    "alloc_MB": 1234.56,
    "gc_time_sec": 0.12,
    "julia_version": "1.11.0",
    "julia_threads": 8,
    "julia_nprocs": 5,
    "julia_nworkers": 4,
    "dionysos_nparts": 4,
    "hostname": "node01",
    "total_memory_GB": 64.0,
    "worker_info": { "2": {"pid": 1234, "host": "node01", "threads": 8}, ... },
    "SLURM_JOB_ID": "12345",
    ...
}
```

## How Parallelism Works in Dionysos

The parameter flow through the Dionysos source:

1. **Environment variables** → `USE_DISTRIBUTED` / `USE_THREADED` in the example script
2. → `MOI.set(optimizer, "distributed"/"threaded", ...)` on `UniformGridAbstraction.Optimizer`
3. → Delegated to `OptimizerEmptyProblem` fields
4. → `MOI.optimize!` reads the fields and calls `compute_abstract_system_from_concrete_system!`
5. → Dispatcher branches based on the two booleans:
   - **Serial**: sequential double loop over (state × input)
   - **Threaded**: `Threads.@threads` parallel loop, results merged per-thread
   - **Distributed**: `Distributed.pmap` partitions the state space across workers
   - **Hybrid**: `pmap` across workers, each worker uses `Threads.@threads`

## Tuning SLURM Parameters

| Parameter | Controls | Maps to Julia |
|-----------|----------|---------------|
| `--nodes` | Physical nodes | Multi-node distribution |
| `--ntasks` | OS processes | `nprocs()` / `nworkers()` |
| `--cpus-per-task` | Cores per process | `Threads.nthreads()` |
| `--mem-per-cpu` | RAM per core | Available to each thread/worker |

**Rule of thumb:** `ntasks × cpus-per-task ≤ total_cores_per_node × nodes`

### Recommended SLURM Settings

| Mode | `--nodes` | `--ntasks` | `--cpus-per-task` | Julia flags |
|------|-----------|------------|-------------------|-------------|
| Serial | 1 | 1 | 1 | (none) |
| Threaded | 1 | 1 | 4–32 | `--threads=$SLURM_CPUS_PER_TASK` |
| Distributed | 1–N | 4–32 | 1 | (workers via `addprocs`) |
| Hybrid | 1–N | 2–8 | 4–8 | `--threads=$SLURM_CPUS_PER_TASK` |
