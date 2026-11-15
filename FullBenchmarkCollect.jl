# Safe to include alongside the development evaluators, or more than once.
if !isdefined(@__MODULE__, :FullBenchmarkRuntime)
    include(joinpath(@__DIR__, "FullBenchmarkRuntime.jl"))
end

prepare_full_benchmark(; kwargs...) = FullBenchmarkRuntime.prepare_full_benchmark(; kwargs...)
run_full_benchmark(; kwargs...) = FullBenchmarkRuntime.run_full_benchmark(; kwargs...)
full_benchmark_status(; kwargs...) = FullBenchmarkRuntime.full_benchmark_status(; kwargs...)

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS == ["--run"]
        run_full_benchmark()
    elseif ARGS == ["--status"]
        println(full_benchmark_status())
    elseif isempty(ARGS) || ARGS == ["--prepare"]
        data = prepare_full_benchmark()
        println("Prepared $(data["prompt_count"]) prompts × $(data["config"]["K"]) queries = $(data["planned_queries"]) planned responses. No API calls made.")
        println("Output: $(FullBenchmarkRuntime.FULL_BENCHMARK_RESULTS_PATH)")
    else
        error("Use --prepare, --run, or --status")
    end
end
