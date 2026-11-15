if !isdefined(@__MODULE__, :VLTLRuntime)
    include(joinpath(@__DIR__, "VLTLRuntime.jl"))
end
prepare_vltl(;kwargs...)=VLTLRuntime.prepare_vltl(;kwargs...)
run_vltl(;kwargs...)=VLTLRuntime.run_vltl(;kwargs...)
vltl_status(;kwargs...)=VLTLRuntime.vltl_status(;kwargs...)
analyze_vltl(;kwargs...)=VLTLRuntime.analyze_vltl(;kwargs...)

if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--run"]
        println(run_vltl())
        analyze_vltl()
    elseif ARGS==["--analyze"]
        analyze_vltl()
    elseif ARGS==["--status"]
        println(vltl_status())
    elseif isempty(ARGS) || ARGS==["--prepare"]
        data=prepare_vltl()
        println("Prepared $(data["prompt_count"]) VLTL prompts, $(data["planned_queries"]) planned responses. No API calls.")
    else
        error("Use --prepare, --run, --analyze, or --status")
    end
end
