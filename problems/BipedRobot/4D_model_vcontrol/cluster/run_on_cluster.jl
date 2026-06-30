# ==============================================================================
# SLURM-array chunk computation for robot abstraction
# ==============================================================================
 
import Pkg
Pkg.instantiate()
 
using Printf
 
global begin
    using Dionysos
    using MathematicalSystems
    using StaticArrays
    using LinearAlgebra
    using JuMP
    using MathOptInterface
    using JLD2
end
 
global const MOI = MathOptInterface
global const DI = Dionysos
global const UT = DI.Utils
global const ST = DI.System
global const MP = DI.Mapping
global const SY = DI.Symbolic
global const OP = DI.Optim
global const AB = OP.Abstraction
 
# ------------------------------------------------------------------------------
# Optimizer builder
# ------------------------------------------------------------------------------
 
function build_optimizer(
    concrete_problem,
    execution_backend,
    state_grid,
    input_grid;
    print_level::Int = 2,
    progress_update_interval::Int = Int(1e2),
)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
 
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
 
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), false)
    # MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), concrete_system.X)
 
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    ) # CENTER_SIMULATION RANDOM_SIMULATION
 
    MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), execution_backend)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("progress_update_interval"),
        progress_update_interval,
    )
    return optimizer
end
 
# ------------------------------------------------------------------------------
# Problem
# ------------------------------------------------------------------------------
 
global begin
    include(joinpath(@__DIR__, "4D_model_vcontrol/cluster/robot_cluster_simulation.jl"))
end
 
# ------------------------------------------------------------------------------
# Chunk parameters
# ------------------------------------------------------------------------------
 
nchunks = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))
chunk_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
 
execution_backend = SY.SlurmArrayBackend(
    nchunks,
    chunk_id,
    outdir,
    :contiguous,
    true, # write_only: do not add transitions to local abstract system
)
 
@info(
    "Starting SLURM abstraction chunk",
    chunk_id = chunk_id,
    nchunks = nchunks,
    outdir = outdir,
)
 
# ------------------------------------------------------------------------------
# Optimizer
# ------------------------------------------------------------------------------
 
global optimizer = build_optimizer(
    concrete_problem,
    execution_backend,
    state_grid,
    input_grid;
    print_level = 2,
    progress_update_interval = Int(1e3),
)
 
# ------------------------------------------------------------------------------
# Build and compute chunk
# ------------------------------------------------------------------------------
 
@printf("Optimizer build time:         %.3f s\n", t_build_optimizer)
 
t_opt_wall = @elapsed MOI.optimize!(optimizer)
 
@printf("Chunk optimize! wall time:    %.3f s\n", t_opt_wall)
 
@info(
    "Finished SLURM abstraction chunk",
    chunk_id = chunk_id,
    nchunks = nchunks,
    elapsed = t_opt_wall,
)