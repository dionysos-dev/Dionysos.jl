include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "non_linear.jl"))
const NL = NonLinear

include(joinpath(@__DIR__, "helpers_nonlinear_bm_vs_mppi.jl"))

function main(cfg::NonlinearBMvsMPPIConfig = NonlinearBMvsMPPIConfig())
    result = run_all_horizons(cfg)
    println("benchmark = nonlinear_bm_vs_mppi")
    println("output_root = ", result.paths.root)
    println("summary_csv = ", joinpath(result.paths.summary_dir, "summary.csv"))
    println("summary_md = ", joinpath(result.paths.summary_dir, "summary.md"))
    println("note = ", joinpath(result.paths.summary_dir, "benchmark_notes.md"))
    println("bm_internal_model = benchmark-local PWA/hybrid surrogate")
    println("mppi_internal_model = true nonlinear discrete rollout")
    println("fairness = same true dynamics/sets/horizon/evaluation/certifier")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
