# ==============================================================================
#  Assemble partition results into a complete abstract system.
#
#  Usage:
#      julia --project assemble_partitions.jl <setup_script> <partitions_dir> <nparts> \
#            [output_dir] [--strategy=roundrobin] [--solve]
#
#  Example (local):
#      julia --project=. parallel_tests/src/assemble_partitions.jl \
#            problems/BipedRobot/6D_model/robot_example_setup.jl \
#            /tmp/partitions 300 /tmp/assembled
#
#  Example (SLURM):
#      See slurm/run_assemble.sh
# ==============================================================================

using Dates
using JLD2
using JSON

# ---------------------------------------------------------------------------
#  Argument parsing
# ---------------------------------------------------------------------------
SOLVE_PROBLEM = false
STRATEGY = :roundrobin
positional = String[]

for a in ARGS
    if a == "--solve"
        global SOLVE_PROBLEM = true
    elseif startswith(a, "--strategy=")
        global STRATEGY = Symbol(match(r"--strategy=(\w+)", a)[1])
    else
        push!(positional, a)
    end
end

if length(positional) < 3
    error(
        "Usage: assemble_partitions.jl <setup_script> <partitions_dir> <nparts> " *
        "[output_dir] [--strategy=roundrobin] [--solve]",
    )
end

const SETUP_SCRIPT = positional[1]
const PARTITIONS_DIR = positional[2]
const NPARTS = parse(Int, positional[3])
const OUTPUT_DIR = length(positional) >= 4 ? positional[4] : PARTITIONS_DIR

mkpath(OUTPUT_DIR)

println("=" ^ 70)
println("  Dionysos Partition Assembler")
println("=" ^ 70)
println("  Setup script  : $(SETUP_SCRIPT)")
println("  Partitions dir: $(PARTITIONS_DIR)")
println("  N parts       : $(NPARTS)")
println("  Output dir    : $(OUTPUT_DIR)")
println("  Solve         : $(SOLVE_PROBLEM)")
println("  Strategy      : $(STRATEGY)")
println("=" ^ 70)
println()

# ---------------------------------------------------------------------------
#  Setup — builds the same empty abstract system
# ---------------------------------------------------------------------------
t_wall_start = time()

println("Loading setup...")
include(SETUP_SCRIPT)
t0 = time()
env = RobotExampleSetup.setup(; partition_strategy = STRATEGY)
t_setup = time() - t0
println("Setup time: $(round(t_setup; digits = 2)) s")

# Convenience aliases
const DI = RobotExampleSetup.DI
const SY = RobotExampleSetup.SY
const UT = RobotExampleSetup.UT
const ST = RobotExampleSetup.ST
const MP = RobotExampleSetup.MP
const OP = RobotExampleSetup.OP
const AB = RobotExampleSetup.AB

# ---------------------------------------------------------------------------
#  Load and merge all partition results
# ---------------------------------------------------------------------------
println("\nLoading partition results...")
all_transitions = Tuple{Int, Int, Int}[]
total_sources = 0
total_compute_sec = 0.0
missing_parts = Int[]
partition_meta = Dict{Int, Dict{String, Any}}()

for i in 1:NPARTS
    path = joinpath(PARTITIONS_DIR, "partition_$(i).jld2")
    if !isfile(path)
        push!(missing_parts, i)
        continue
    end

    data = load(path)
    transitions = data["transitions"]::Vector{Tuple{Int, Int, Int}}
    meta = data["metadata"]::Dict{String, Any}

    append!(all_transitions, transitions)
    total_sources += meta["n_source_states"]
    total_compute_sec += meta["elapsed_compute_sec"]
    partition_meta[i] = meta

    println(
        "  Partition $(lpad(i, ndigits(NPARTS))): " *
        "$(lpad(meta["n_transitions"], 8)) transitions, " *
        "$(lpad(meta["n_source_states"], 6)) states, " *
        "$(round(meta["elapsed_compute_sec"]; digits = 1))s",
    )
end

if !isempty(missing_parts)
    @warn "Missing partition files" count = length(missing_parts) indices = missing_parts
end

println("\nTotal transitions loaded : $(length(all_transitions))")
println("Total source states      : $(total_sources)")
println("Sum of compute times     : $(round(total_compute_sec; digits = 1)) s")

# ---------------------------------------------------------------------------
#  Add transitions to the abstract system
# ---------------------------------------------------------------------------
println("\nAdding transitions to abstract system...")
t0 = time()
SY.add_transitions!(env.abstract_system, all_transitions)
elapsed_add = time() - t0
println("add_transitions! time: $(round(elapsed_add; digits = 2)) s")
println("Total transitions in automaton: $(SY.get_n_transitions(env.abstract_system))")

# ---------------------------------------------------------------------------
#  Save assembled abstraction via the optimizer's JLD2 export
# ---------------------------------------------------------------------------
# Build an Optimizer wrapper so we can use the standard export_abstraction_jld2
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
# Create a minimal abstraction_solver and inject the abstract system
if optimizer.abstraction_solver === nothing
    optimizer.abstraction_solver = AB.UniformGridAbstraction.OptimizerEmptyProblem()
end
optimizer.abstraction_solver.abstract_system = env.abstract_system

assembled_path = joinpath(OUTPUT_DIR, "assembled_abstraction.jld2")
AB.UniformGridAbstraction.export_abstraction_jld2(optimizer, assembled_path)
println("Assembled abstraction saved to: $(assembled_path)")

# ---------------------------------------------------------------------------
#  Save assembly metadata as JSON
# ---------------------------------------------------------------------------
assembly_meta = Dict{String, Any}(
    "nparts" => NPARTS,
    "n_loaded" => NPARTS - length(missing_parts),
    "n_missing" => length(missing_parts),
    "missing_parts" => missing_parts,
    "total_transitions" => length(all_transitions),
    "total_source_states" => total_sources,
    "sum_compute_sec" => total_compute_sec,
    "add_transitions_sec" => elapsed_add,
    "setup_time_sec" => t_setup,
    "total_wall_clock_sec" => time() - t_wall_start,
    "timestamp" => string(Dates.now()),
    "hostname" => gethostname(),
    "strategy" => string(STRATEGY),
)

assembly_meta_path = joinpath(OUTPUT_DIR, "assembly_metadata.json")
open(assembly_meta_path, "w") do io
    return JSON.print(io, assembly_meta, 4)
end
println("Assembly metadata saved to: $(assembly_meta_path)")

# ---------------------------------------------------------------------------
#  Optionally solve the abstract optimal-control problem
# ---------------------------------------------------------------------------
if SOLVE_PROBLEM
    using MathematicalSystems
    import StaticArrays: SVector

    println("\n" * "=" ^ 70)
    println("  Solving abstract optimal-control problem")
    println("=" ^ 70)

    n_state = env.n_state
    concrete_system = env.concrete_system

    # --- First step ---
    x0 = SVector{n_state, Float64}(zeros(n_state))
    t_low = SVector{n_state, Float64}([-12π / 180, 7π / 180, 8π / 180, -0.75, -0.30, -0.30])
    t_high = SVector{n_state, Float64}([-8π / 180, 9π / 180, 12π / 180, 0.30, 0.75, 0.75])

    I = UT.HyperRectangle(x0, x0)
    T = UT.HyperRectangle(t_low, t_high)

    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        I,
        T,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )

    # Set the control problem on the optimizer (abstraction is already cached)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

    MOI.optimize!(optimizer)

    t_solve = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
    println("Solve time: $(t_solve) s")

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    stopfun = x -> (x ∈ problem.target_set)

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        concrete_system,
        controller,
        x0,
        300;
        stopping = stopfun,
    )

    traj_path = joinpath(OUTPUT_DIR, "trajectory_step1_assembled.jld2")
    @save traj_path x_traj u_traj
    println("Trajectory saved to: $(traj_path)")
    println("Final state: ", x_traj.seq[end])
end

# ---------------------------------------------------------------------------
total_wall = time() - t_wall_start
println("\n" * "=" ^ 70)
println("  Assembly complete!")
println("  Total wall-clock: $(round(total_wall; digits = 2)) s")
println("=" ^ 70)
