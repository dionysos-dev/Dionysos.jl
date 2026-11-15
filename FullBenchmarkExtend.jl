# Adding the extension method to the already-loaded runtime also works in an
# existing REPL; the shared module and its constants do not need to be reloaded.
include(joinpath(@__DIR__, "FullBenchmarkCollect.jl"))

if !isdefined(FullBenchmarkRuntime, :extend_full_benchmark)
@eval FullBenchmarkRuntime begin
    """Seed a larger-R JSON from a completed smaller-R run, without API calls.

    `K` is the collector's number of independent samples (R), not the number of
    semantic clusters. The original JSON remains unchanged. Repeating this
    operation validates and returns an already-seeded target, including partial
    progress from queries beyond the original R.
    """
    function extend_full_benchmark(;
        source_path=joinpath(@__DIR__, "Results", "deepseek-v4-pro-task-specific-full-T0.7-K7.json"),
        results_path=joinpath(@__DIR__, "Results", "deepseek-v4-pro-task-specific-full-T0.7-K10.json"),
        K=10, dataset_path=DEFAULT_DATASET_PATH,
        development_manifest_path=joinpath(@__DIR__, "Results", "task-specific-temperature-prompt-v2", "manifest.json"))
        K isa Integer && !(K isa Bool) && K>0 || throw(ArgumentError("K must be a positive integer"))
        source=abspath(source_path)
        target=abspath(results_path)
        isfile(source) || throw(ArgumentError("Source checkpoint not found: $source"))
        (source==target || (isfile(target) && realpath(source)==realpath(target))) &&
            throw(ArgumentError("Use a separate target JSON to preserve the original experiment"))
        # Lock both paths in a fixed order, including the read-only source, so a
        # collector or another extension cannot change a checkpoint during copying.
        paths=sort([source,target])
        return with_collection_lock(paths[1]) do
            with_collection_lock(paths[2]) do
                original=mutable_json(JSON3.read(read(source,String)))
                config=original["config"]
                old_K=config["K"]
                old_K isa Integer && !(old_K isa Bool) && K>old_K ||
                    throw(ArgumentError("New K=$K must exceed the source sample count R=$old_K"))
                settings=(dataset_path=String(dataset_path), K=Int(old_K),
                    temperature=Float64(config["temperature"]), model=String(config["model"]),
                    input_fields=String.(config["input_fields"]), max_tokens=Int(config["max_tokens"]),
                    api_url=String(config["api_url"]), api_key_env=String(config["api_key_env"]),
                    development_manifest_path=development_manifest_path)
                # Reuse the normal resume validator for model, dataset, messages,
                # allowed APs, token settings, metadata, and every query index.
                original=prepare_collection(source; settings...)
                collection_summary(original)["complete"] ||
                    error("Finish the source R=$old_K collection before extending it")
                target_settings=merge(settings,(K=Int(K),))
                if isfile(target)
                    existing=prepare_collection(target; target_settings...)
                    for (old,p) in zip(original["prompts"],existing["prompts"])
                        length(p["samples"])>=old_K && p["samples"][1:old_K]==old["samples"] ||
                            error("Target prompt $(p["prompt_index"]) does not preserve the source responses; use a new results_path")
                    end
                    println("Using existing extension: $(collection_summary(existing)["remaining_queries"]) responses remain. No API calls made.")
                    return existing
                end
                extended=deepcopy(original)
                extended["config"]["K"]=Int(K)
                extended["config"]["repetitions"]=Int(K)
                extended["planned_queries"]=extended["prompt_count"]*K
                extended["dataset_path"]=abspath(dataset_path)
                extended["created_at_utc"]=utc_timestamp()
                extended["run_state"]="prepared"
                history=get!(extended,"extension_history",Any[])
                push!(history,OrderedDict{String,Any}(
                    "source_path"=>source,"source_sha256"=>bytes2hex(sha256(read(source))),
                    "source_repetitions"=>old_K,"target_repetitions"=>K,
                    "inherited_responses"=>collection_summary(original)["saved_responses"],
                    "source_created_at_utc"=>original["created_at_utc"],"extended_at_utc"=>utc_timestamp()))
                summary=checkpoint!(extended,target)
                println("Preserved $(summary["saved_responses"]) responses; $(summary["remaining_queries"]) additional responses needed for R=$K. No API calls made.")
                println("Target: $target")
                return extended
            end
        end
    end
end
end # runtime method guard

extend_full_benchmark(;kwargs...)=FullBenchmarkRuntime.extend_full_benchmark(;kwargs...)

if abspath(PROGRAM_FILE)==@__FILE__
    isempty(ARGS) || ARGS==["--prepare"] || error("Use --prepare (copies saved responses without API calls)")
    extend_full_benchmark()
end
