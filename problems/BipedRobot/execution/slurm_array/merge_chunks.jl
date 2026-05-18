include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())

nchunks = parse(Int, get(ENV, "DIONYSOS_NCHUNKS", "1"))
simplify = parse(Float64, get(ENV, "DIONYSOS_SIMPLIFY", "3.0"))
tstep = parse(Float64, get(ENV, "DIONYSOS_TSTEP", "0.1"))

@info(
    "Preparing empty abstraction for merge",
    model = selected_robot_model(),
    nchunks = nchunks,
    outdir = outdir,
    outfile = outfile,
    simplify = simplify,
    tstep = tstep,
)

execution_backend = SY.SlurmArrayBackend(nchunks, 1, outdir, :contiguous, true)

optimizer = build_robot_abstraction_optimizer(;
    execution_backend = execution_backend,
    simplify = simplify,
    tstep = tstep,
    print_level = 1,
)

abstract_system = build_empty_abstraction_for_optimizer!(optimizer)

@info "Merging transition chunks" nchunks outdir

merge_time = @elapsed SY.merge_transition_chunks!(
    abstract_system,
    outdir;
    nchunks = nchunks,
    print_level = 1,
)

@info "Merge finished" elapsed = merge_time

save_optimizer(outfile, optimizer)

@info "Saved merged optimizer" outfile
