# ==============================================================================
#  Compute transitions for a single partition of the abstract state space.
#
#  Usage:
#      julia --project partition_runner.jl <setup_script> <partition_idx> <nparts> <output_dir> \
#            [--strategy=roundrobin]
#
#  Example (local):
#      julia --project=. parallel_tests/src/partition_runner.jl \
#            problems/BipedRobot/6D_model/robot_example_setup.jl \
#            1 300 /tmp/partitions
#
#  Example (SLURM array):
#      See slurm/run_array_partitions.sh
# ==============================================================================

using Dates
using JLD2

# ---------------------------------------------------------------------------
#  Argument parsing
# ---------------------------------------------------------------------------
if length(ARGS) < 4
    error(
        "Usage: partition_runner.jl <setup_script> <partition_idx> <nparts> <output_dir> " *
        "[--strategy=roundrobin]",
    )
end

const SETUP_SCRIPT  = ARGS[1]
const PARTITION_IDX = parse(Int, ARGS[2])
const NPARTS        = parse(Int, ARGS[3])
const OUTPUT_DIR    = ARGS[4]

STRATEGY = :roundrobin
for a in ARGS[5:end]
    m = match(r"--strategy=(\w+)", a)
    m !== nothing && (global STRATEGY = Symbol(m[1]))
end

mkpath(OUTPUT_DIR)

println("=" ^ 70)
println("  Dionysos Partition Runner")
println("=" ^ 70)
println("  Setup script : $(SETUP_SCRIPT)")
println("  Partition    : $(PARTITION_IDX) / $(NPARTS)")
println("  Strategy     : $(STRATEGY)")
println("  Output dir   : $(OUTPUT_DIR)")
println("  PID          : $(getpid())")
println("  Host         : $(gethostname())")
println("  Threads      : $(Threads.nthreads())")
println("  Julia        : $(VERSION)")
println("=" ^ 70)
println()

# ---------------------------------------------------------------------------
#  Setup — loads packages, builds abstract system (no transitions)
# ---------------------------------------------------------------------------
t_wall_start = time()

println("Loading setup...")
include(SETUP_SCRIPT)
t0 = time()
env = RobotExampleSetup.setup(; partition_strategy = STRATEGY)
t_setup = time() - t0
println("Setup time: $(round(t_setup; digits = 2)) s")

# Convenience aliases
const SY = RobotExampleSetup.SY
const MP = RobotExampleSetup.MP

# ---------------------------------------------------------------------------
#  Partition the source states and select this partition
# ---------------------------------------------------------------------------
println("\nPartitioning source states ($(NPARTS) parts, $(STRATEGY))...")
t0 = time()
parts = SY.partition_source_states(
    env.abstract_system, NPARTS; strategy = STRATEGY,
)
t_partition = time() - t0

n_local_states = length(collect(
    MP.enum_states(parts[PARTITION_IDX], SY.get_state_mapping(env.abstract_system)),
))
println("Partition time : $(round(t_partition; digits = 2)) s")
println("Local states   : $(n_local_states)")

# Build a LocalGridBasedSymbolicModel for this partition
local_symmodel = SY.LocalGridBasedSymbolicModel(
    env.abstract_system, parts[PARTITION_IDX],
)

# ---------------------------------------------------------------------------
#  Compute transitions
# ---------------------------------------------------------------------------
println("\nComputing transitions for partition $(PARTITION_IDX)...")
RobotExampleSetup.RobotProblem.reset_rbd_timing!()
t0 = time()
transitions = SY.collect_abstract_transitions(
    local_symmodel,
    env.system_approximation;
    print_level      = 0,
    state_filter      = env.state_filter,
    state_input_filter = env.state_input_filter,
)
elapsed_compute = time() - t0

rbd_stats = RobotExampleSetup.RobotProblem.get_rbd_timing()

println("Compute time   : $(round(elapsed_compute; digits = 2)) s")
println("Transitions    : $(length(transitions))")
println("RBD sim time   : $(round(rbd_stats.time_sec; digits = 2)) s")
println("RBD calls      : $(rbd_stats.call_count)")

# ---------------------------------------------------------------------------
#  Save results to JLD2
# ---------------------------------------------------------------------------
output_file = joinpath(OUTPUT_DIR, "partition_$(PARTITION_IDX).jld2")

metadata = Dict{String, Any}(
    "partition_idx"       => PARTITION_IDX,
    "nparts"              => NPARTS,
    "strategy"            => string(STRATEGY),
    "n_source_states"     => n_local_states,
    "n_transitions"       => length(transitions),
    "elapsed_compute_sec" => elapsed_compute,
    "setup_time_sec"      => t_setup,
    "partition_time_sec"  => t_partition,
    "rbd_time_sec"        => rbd_stats.time_sec,
    "rbd_call_count"      => rbd_stats.call_count,
    "hostname"            => gethostname(),
    "pid"                 => getpid(),
    "timestamp"           => string(Dates.now()),
    "julia_version"       => string(VERSION),
    "julia_threads"       => Threads.nthreads(),
)

@save output_file transitions metadata
println("\nSaved: $(output_file)")

total_wall = time() - t_wall_start
println("\nTotal wall-clock: $(round(total_wall; digits = 2)) s")
