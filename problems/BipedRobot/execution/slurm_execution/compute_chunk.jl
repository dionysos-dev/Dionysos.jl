include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
nchunks = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))
chunk_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

simplify = parse(Float64, get(ENV, "DIONYSOS_SIMPLIFY", "3.0"))
tstep = parse(Float64, get(ENV, "DIONYSOS_TSTEP", "0.1"))

@info(
    "Starting SLURM abstraction chunk",
    model = selected_robot_model(),
    chunk_id = chunk_id,
    nchunks = nchunks,
    outdir = outdir,
    simplify = simplify,
    tstep = tstep,
)

execution_backend = SY.SlurmArrayBackend(nchunks, chunk_id, outdir, :contiguous, true)

optimizer = build_robot_abstraction_optimizer(;
    execution_backend = execution_backend,
    simplify = simplify,
    tstep = tstep,
    print_level = 2,
)

elapsed = @elapsed MOI.optimize!(optimizer)

@info(
    "Finished SLURM abstraction chunk",
    chunk_id = chunk_id,
    nchunks = nchunks,
    elapsed = elapsed,
)
