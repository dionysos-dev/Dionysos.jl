#!/usr/bin/env julia
"""
Wrapper that runs a Dionysos example script, captures wall-clock time and
system/environment metrics, and writes a JSON results file.

Usage:
    julia --project=<PROJECT> run_and_collect.jl <example.jl>

The example path must be absolute or relative to the current directory.
Environment variables DIONYSOS_DISTRIBUTED, DIONYSOS_THREADED, DIONYSOS_NPARTS,
and DIONYSOS_OUTDIR are read and recorded in the metrics.
"""

using JSON
using Dates
using Distributed

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
example_path = ARGS[1]
if !isfile(example_path)
    error("Example script not found: $example_path")
end

outdir = get(ENV, "DIONYSOS_OUTDIR", dirname(example_path))
mkpath(outdir)

is_distributed = lowercase(get(ENV, "DIONYSOS_DISTRIBUTED", "false")) == "true"
is_threaded    = lowercase(get(ENV, "DIONYSOS_THREADED", "false")) == "true"
nparts         = parse(Int, get(ENV, "DIONYSOS_NPARTS", "1"))

mode = if is_distributed && is_threaded
    "hybrid"
elseif is_distributed
    "distributed"
elseif is_threaded
    "threaded"
else
    "serial"
end

println("=" ^ 60)
println("  Dionysos Benchmark Wrapper")
println("=" ^ 60)
println("  Example:      $(example_path)")
println("  Mode:         $(mode)")
println("  Distributed:  $(is_distributed)")
println("  Threaded:     $(is_threaded)")
println("  NPARTS:       $(nparts)")
println("  Threads:      $(Threads.nthreads())")
println("  PID:          $(getpid())")
println("  Host:         $(gethostname())")
println("  Output:       $(outdir)")
println("=" ^ 60)
println()

# ---------------------------------------------------------------------------
# Run the example and measure time
# ---------------------------------------------------------------------------
t0 = time()
stats = @timed include(example_path)
elapsed = stats.time
alloc_bytes = stats.bytes
gc_time = stats.gctime
total_time = time() - t0

println()
println("=" ^ 60)
println("  Benchmark finished")
println("  Wall-clock:   $(round(total_time; digits=4)) s")
println("  @timed:       $(round(elapsed; digits=4)) s")
println("  Alloc:        $(round(alloc_bytes / 1e6; digits=2)) MB")
println("  GC time:      $(round(gc_time; digits=4)) s")
println("=" ^ 60)

# ---------------------------------------------------------------------------
# Collect worker info (if distributed)
# ---------------------------------------------------------------------------
worker_info = Dict{String, Any}()
if nworkers() > 1
    for wid in workers()
        info = @fetchfrom wid begin
            Dict{String, Any}(
                "pid"     => getpid(),
                "host"    => gethostname(),
                "threads" => Threads.nthreads(),
            )
        end
        worker_info[string(wid)] = info
    end
end

# ---------------------------------------------------------------------------
# Build metrics
# ---------------------------------------------------------------------------
metrics = Dict{String, Any}(
    "mode"              => mode,
    "example"           => basename(example_path),
    "example_path"      => example_path,
    "timestamp"         => string(Dates.now()),
    "hostname"          => gethostname(),
    "pid"               => getpid(),
    "elapsed_sec"       => round(elapsed; digits=6),
    "wall_clock_sec"    => round(total_time; digits=6),
    "alloc_MB"          => round(alloc_bytes / 1e6; digits=2),
    "gc_time_sec"       => round(gc_time; digits=6),
    "julia_version"     => string(VERSION),
    "julia_threads"     => Threads.nthreads(),
    "julia_nprocs"      => nprocs(),
    "julia_nworkers"    => nworkers(),
    "dionysos_nparts"   => nparts,
    "total_memory_GB"   => round(Sys.total_memory() / 1e9; digits=2),
    "free_memory_GB"    => round(Sys.free_memory() / 1e9; digits=2),
)

if !isempty(worker_info)
    metrics["worker_info"] = worker_info
    hosts = unique([info["host"] for (_, info) in worker_info])
    metrics["nodes_used"] = hosts
end

# Capture SLURM environment variables
for var in ["SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_NNODES", "SLURM_NTASKS",
            "SLURM_CPUS_PER_TASK", "SLURM_NODELIST", "SLURM_MEM_PER_CPU"]
    val = get(ENV, var, nothing)
    if val !== nothing
        metrics[var] = val
    end
end

# ---------------------------------------------------------------------------
# Write JSON
# ---------------------------------------------------------------------------
run_id = get(ENV, "SLURM_JOB_ID", string(getpid()))
json_path = joinpath(outdir, "metrics_$(mode)_$(run_id).json")
mkpath(dirname(json_path))
open(json_path, "w") do io
    JSON.print(io, metrics, 4)
end
println("Metrics saved to: $(json_path)")
