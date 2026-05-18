import Pkg
Pkg.instantiate()

using Distributed
using Printf

const NWORKERS = parse(Int, get(ENV, "DIONYSOS_NWORKERS", "4"))
const USE_THREADED_PER_WORKER =
    lowercase(get(ENV, "DIONYSOS_THREADED_PER_WORKER", "false")) == "true"

const DISTRIBUTED_NPARTS = parse(Int, get(ENV, "DIONYSOS_NPARTS", string(NWORKERS)))
const PARTITION_STRATEGY = Symbol(get(ENV, "DIONYSOS_PARTITION_STRATEGY", "contiguous"))

const SIMPLIFY = parse(Float64, get(ENV, "DIONYSOS_SIMPLIFY", "3.0"))
const TSTEP = parse(Float64, get(ENV, "DIONYSOS_TSTEP", "0.1"))

const BIPED_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const PROJECT_DIR = BIPED_ROOT

if length(workers()) < NWORKERS
    addprocs(NWORKERS - length(workers()) + 1; exeflags = `--project=$(PROJECT_DIR)`)
end

@everywhere using Dionysos

include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

robot_problem_path = selected_robot_problem_path()

@everywhere include($robot_problem_path)
@everywhere using .RobotProblem

@info(
    "Starting Julia distributed abstraction",
    model = selected_robot_model(),
    nworkers = length(workers()),
    nparts = DISTRIBUTED_NPARTS,
    partition_strategy = PARTITION_STRATEGY,
    threaded_per_worker = USE_THREADED_PER_WORKER,
    simplify = SIMPLIFY,
    tstep = TSTEP,
)

execution_backend = SY.JuliaDistributedBackend(
    nothing,
    DISTRIBUTED_NPARTS,
    PARTITION_STRATEGY,
    USE_THREADED_PER_WORKER,
    true,
)

optimizer = build_robot_abstraction_optimizer(;
    execution_backend = execution_backend,
    simplify = SIMPLIFY,
    tstep = TSTEP,
    print_level = 2,
)

elapsed = @elapsed MOI.optimize!(optimizer)

construct_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))

@printf("Abstraction optimize! wall time: %.3f s\n", elapsed)
@printf("Abstraction reported construction time: %.3f s\n", construct_time)

outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())
save_optimizer(outfile, optimizer)

@info "Saved abstraction" outfile
