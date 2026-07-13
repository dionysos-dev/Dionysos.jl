# ==============================================================================
# Merge SLURM-array transition chunks into one abstraction
# ==============================================================================
 
import Pkg
Pkg.instantiate()
 
using Printf
 
using Dionysos
using MathematicalSystems
using StaticArrays
using LinearAlgebra
using JuMP
using MathOptInterface
using JLD2
 
 
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

    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.1)
 
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

function build_empty_abstraction_for_optimizer!(optimizer)
    abstraction_solver = optimizer.abstraction_solver
    abstraction_solver === nothing && error("Optimizer has no abstraction_solver.")
 
    return AB.UniformGridAbstraction.build_empty_abstraction!(abstraction_solver)
end

function save_optimizer(filename::AbstractString, optimizer)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        return file["optimizer"] = optimizer
    end
    return filename
end

# ------------------------------------------------------------------------------
# Problem
# ------------------------------------------------------------------------------

include(joinpath(@__DIR__, "robot_create_problem.jl"))

concrete_problem = DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)
 
# ------------------------------------------------------------------------------
# Slurm Parameters
# ------------------------------------------------------------------------------
 
nchunks = 100
outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", "")
outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", "")

execution_backend = SY.SlurmArrayBackend(
    nchunks,
    1,
    outdir,
    :contiguous,
    true
)
 
# ------------------------------------------------------------------------------
# Optimizer (same as the one for computing the transition chunks)
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
# Merge
# ------------------------------------------------------------------------------

@info(
    "Preparing empty abstraction for merge",
    nchunks = nchunks,
    outdir = outdir,
    outfile = outfile,
)
 
abstract_system = build_empty_abstraction_for_optimizer!(optimizer)

@info "Merging transition chunks" nchunks outdir
 
t_merge = @elapsed SY.merge_transition_chunks!(
    abstract_system,
    outdir;
    nchunks = nchunks,
    print_level = 1,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
 
@printf("Merge time:                   %.3f s\n", t_merge)
 
t_save = @elapsed JLD2.jldsave(outfile; optimizer = optimizer) #save_optimizer(outfile, optimizer)
 
@printf("Save time:                    %.3f s\n", t_save)
 
@info "Saved merged optimizer" outfile