module FullBenchmarkRuntime

# Isolate the shared evaluator from constants already loaded in the user's REPL.
# Collection does not call its parsing, Spot, clustering, or synthesis functions.
include(joinpath(@__DIR__, "taskSpescific_evaluate.jl"))
using SHA

const FULL_BENCHMARK_RESULTS_PATH = joinpath(@__DIR__, "Results",
    "deepseek-v4-pro-task-specific-full-T0.7-K7.json")

utc_timestamp() = string(now(UTC)) * "Z"
mutable_json(x::AbstractDict) = OrderedDict{String,Any}(String(k) => mutable_json(v) for (k,v) in pairs(x))
mutable_json(x::AbstractVector) = Any[mutable_json(v) for v in x]
mutable_json(x) = x

function with_collection_lock(f, results_path)
    ensure_parent_directory(results_path)
    lock_path = results_path * ".lock"
    try
        mkdir(lock_path)
    catch
        isdir(lock_path) || rethrow()
        error("Collection is locked: $lock_path. Do not run two collectors against the same JSON. " *
              "If a previous Julia process was killed, remove this directory only after verifying it has stopped.")
    end
    try
        write(joinpath(lock_path, "owner.txt"), "PID=$(getpid())\nStarted=$(utc_timestamp())\n")
        return f()
    finally
        rm(lock_path; recursive=true)
    end
end

function collection_plan(; dataset_path=DEFAULT_DATASET_PATH, K=7, temperature=0.7,
    model=DEFAULT_MODEL, input_fields=DEFAULT_INPUT_FIELDS, max_tokens=DEFAULT_MAX_TOKENS,
    api_url=DEFAULT_API_URL, api_key_env=DEFAULT_API_KEY_ENV,
    development_manifest_path=joinpath(@__DIR__, "Results", "task-specific-temperature-prompt-v2", "manifest.json"))
    K isa Integer && !(K isa Bool) && K > 0 || throw(ArgumentError("K must be a positive integer"))
    isfinite(temperature) && 0 <= temperature <= 2 || throw(ArgumentError("Invalid temperature"))
    max_tokens isa Integer && max_tokens > 0 || throw(ArgumentError("max_tokens must be positive"))
    fields = String.(collect(input_fields))
    !isempty(fields) && length(unique(fields)) == length(fields) ||
        throw(ArgumentError("input_fields must be nonempty and unique"))
    records = load_dataset(String(dataset_path))
    isempty(records) && throw(ArgumentError("Dataset is empty"))
    development = Dict{Tuple{String,String},Int}()
    if development_manifest_path !== nothing && isfile(development_manifest_path)
        dev = JSON3.read(read(development_manifest_path, String))
        for p in dev["prompts"]
            development[(string(p["record_id"]), String(p["input_field"]))] = Int(p["prompt_index"])
        end
    end
    config = OrderedDict{String,Any}(
        "dataset_sha256" => bytes2hex(sha256(read(dataset_path))),
        "provider" => "deepseek", "model" => String(model),
        "temperature" => Float64(temperature), "K" => Int(K), "repetitions" => Int(K),
        "input_fields" => fields, "max_tokens" => Int(max_tokens),
        "api_url" => String(api_url), "api_key_env" => String(api_key_env),
        "thinking_type" => "disabled", "prompt_version" => TASK_SPECIFIC_PROMPT_VERSION,
        "provider_adapter_version" => TASK_PROVIDER_ADAPTER_VERSION,
        "system_prompt" => TASK_TRANSLATION_SYSTEM_PROMPT,
        "atomic_propositions" => copy(TASK_SPECIFIC_ATOMIC_PROPOSITIONS),
        "response_policy" => "retain_all_returned_outputs_without_logical_filtering",
    )
    prompts = Any[]
    seen_ids = Set{String}()
    for record in records
        haskey(record, "id") || throw(ArgumentError("A benchmark record has no id"))
        rid = string(record["id"])
        rid in seen_ids && throw(ArgumentError("Duplicate record id: $rid"))
        push!(seen_ids, rid)
        for field in fields
            text = get(record, field, nothing)
            text isa AbstractString && !isempty(strip(text)) ||
                throw(ArgumentError("Record $rid has a missing/empty field: $field"))
            mapping = mutable_json(get(record, "task_specific_ap_mapping", nothing))
            prompt = build_task_specific_translation_prompt(text; ap_mapping=mapping)
            push!(prompts, OrderedDict{String,Any}(
                "prompt_index" => length(prompts)+1, "record_id" => record["id"],
                "template_id" => get(record, "template_id", nothing), "input_field" => field,
                "source_text" => String(text), "ground_truth_ltl" => get(record, "LTL", nothing),
                "task_specific_ap_mapping" => mapping, "translation_prompt" => prompt,
                "development_prompt_index" => get(development, (rid, field), nothing),
                "samples" => Any[], "errors" => Any[],
            ))
        end
    end
    return OrderedDict{String,Any}(
        "schema_version" => 1, "stage" => "LLM_collection_only",
        "dataset_path" => abspath(dataset_path), "config" => config,
        "record_count" => length(records), "prompt_count" => length(prompts),
        "planned_queries" => length(prompts)*K,
        "created_at_utc" => utc_timestamp(), "updated_at_utc" => utc_timestamp(),
        "run_state" => "prepared", "prompts" => prompts,
    )
end

function collection_summary(data)
    K = data["config"]["K"]
    saved = complete = errors = output_tokens = input_tokens = missing_usage = 0
    for p in data["prompts"]
        saved += length(p["samples"])
        complete += length(p["samples"]) == K
        errors += length(p["errors"])
        for s in p["samples"]
            usage = s["usage"]
            ot = usage === nothing ? nothing : get(usage, "completion_tokens", nothing)
            it = usage === nothing ? nothing : get(usage, "prompt_tokens", nothing)
            if ot === nothing
                missing_usage += 1
            else
                output_tokens += Int(ot)
            end
            it === nothing || (input_tokens += Int(it))
        end
    end
    return OrderedDict{String,Any}(
        "prompt_count" => data["prompt_count"], "completed_prompts" => complete,
        "planned_queries" => data["planned_queries"], "saved_responses" => saved,
        "remaining_queries" => data["planned_queries"]-saved,
        "complete" => saved == data["planned_queries"], "failed_attempts" => errors,
        "total_output_tokens" => output_tokens, "total_input_tokens" => input_tokens,
        "responses_missing_output_token_usage" => missing_usage,
        "mean_output_tokens" => saved == missing_usage ? nothing : output_tokens/(saved-missing_usage),
    )
end

function checkpoint!(data, results_path)
    data["updated_at_utc"] = utc_timestamp()
    data["summary"] = collection_summary(data)
    # Atomic replacement: an interrupted write leaves the previous valid JSON intact.
    save_results(data, results_path)
    return data["summary"]
end

function prepare_collection(results_path; kwargs...)
    expected = collection_plan(; kwargs...)
    if !isfile(results_path)
        checkpoint!(expected, results_path)
        return expected
    end
    # A corrupt or incompatible checkpoint is an error, never a reason to start over.
    data = mutable_json(JSON3.read(read(results_path, String)))
    get(data, "schema_version", nothing) == 1 && get(data, "stage", nothing) == "LLM_collection_only" ||
        error("Not a full-benchmark collection checkpoint: $results_path")
    data["config"] == expected["config"] || error(
        "Checkpoint configuration/dataset differs. Resume with the original settings or choose a new results_path.")
    for key in ("record_count", "prompt_count", "planned_queries")
        data[key] == expected[key] || error("Checkpoint $key differs from the benchmark")
    end
    length(data["prompts"]) == length(expected["prompts"]) || error("Checkpoint prompt count differs")
    K = data["config"]["K"]
    for (saved, planned) in zip(data["prompts"], expected["prompts"])
        for key in keys(planned)
            # Development tags are annotations; the source manifest need not travel with the checkpoint.
            key in ("samples", "errors", "development_prompt_index") && continue
            saved[key] == planned[key] || error("Checkpoint prompt $(planned["prompt_index"]) differs at $key")
        end
        samples = saved["samples"]
        length(samples) <= K || error("Too many saved queries for prompt $(saved["prompt_index"])")
        for (j, sample) in enumerate(samples)
            sample["query_index"] == j && sample["repetition"] == j && sample["status"] == "returned" ||
                error("Invalid saved query sequence for prompt $(saved["prompt_index"])")
            sample["formula"] isa AbstractString && sample["raw_output"] isa AbstractString ||
                error("Invalid saved output for prompt $(saved["prompt_index"]), query $j")
            haskey(sample, "usage") || error("Saved response is missing usage metadata")
        end
        saved["errors"] isa AbstractVector || error("Invalid saved errors array")
    end
    return data
end

"""Create or validate the one-file experiment without API calls or logical checks."""
function prepare_full_benchmark(; results_path=FULL_BENCHMARK_RESULTS_PATH, kwargs...)
    path = abspath(results_path)
    return with_collection_lock(path) do
        prepare_collection(path; kwargs...)
    end
end

"""Read progress from the single saved JSON; makes no API calls."""
function full_benchmark_status(; results_path=FULL_BENCHMARK_RESULTS_PATH)
    data = JSON3.read(read(results_path, String))
    return collection_summary(data)
end

"""Collect K independent responses per record/field, resuming saved query indices.

All returned outputs count toward K, including invalid, empty, and truncated ones.
HTTP failures stop the run, are saved separately, and are retried on manual resume.
Checkpoint every `checkpoint_every` responses or `checkpoint_seconds` seconds,
and on normal completion, request failure, or a handled interruption.
"""
function run_full_benchmark(; results_path=FULL_BENCHMARK_RESULTS_PATH, workers=4,
    checkpoint_every=50, checkpoint_seconds=30.0, request_fn=request_translation,
    verbose=true, kwargs...)
    workers isa Integer && workers > 0 || throw(ArgumentError("workers must be positive"))
    checkpoint_every isa Integer && checkpoint_every > 0 || throw(ArgumentError("checkpoint_every must be positive"))
    isfinite(checkpoint_seconds) && checkpoint_seconds > 0 || throw(ArgumentError("checkpoint_seconds must be positive"))
    path = abspath(results_path)
    return with_collection_lock(path) do
        data = prepare_collection(path; kwargs...)
        config = data["config"]
        summary = collection_summary(data)
        summary["complete"] && return summary
        # Mock request functions do not need an API key; production always checks it first.
        request_fn === request_translation && get_api_key(config["api_key_env"])
        pending = [i for (i,p) in enumerate(data["prompts"]) if length(p["samples"]) < config["K"]]
        jobs = Channel{Int}(length(pending))
        for i in pending
            put!(jobs, i)
        end
        close(jobs)
        state_lock = ReentrantLock()
        stopped = Ref(false)
        first_error = Ref{Union{Nothing,String}}(nothing)
        dirty = Ref(0)
        last_save = Ref(time())
        tasks = Task[]
        data["run_state"] = "running"
        checkpoint!(data, path)
        verbose && println("Collecting $(summary["remaining_queries"]) remaining responses; K=$(config["K"]), T=$(config["temperature"]), model=$(config["model"]).")
        try
            @sync for _ in 1:min(workers, length(pending))
                task = @async try
                  for index in jobs
                    stopped[] && break
                    p = data["prompts"][index]
                    for j in (length(p["samples"])+1):config["K"]
                        stopped[] && break
                        started = time()
                        sample = try
                            response = request_fn(p["source_text"]; provider=:deepseek,
                                model=config["model"], api_url=config["api_url"], api_key_env=config["api_key_env"],
                                include_temperature=true, temperature=Float64(config["temperature"]),
                                max_tokens=Int(config["max_tokens"]), deepseek_thinking_type="disabled",
                                return_metadata=true, ap_mapping=p["task_specific_ap_mapping"])
                            OrderedDict{String,Any}(
                                "query_index" => j, "repetition" => j, "status" => "returned",
                                "formula" => String(response.formula), "raw_output" => String(response.raw_output),
                                "usage" => mutable_json(response.usage), "returned_model" => response.returned_model,
                                "finish_reason" => response.finish_reason, "elapsed_seconds" => time()-started,
                                "received_at_utc" => utc_timestamp(),
                            )
                        catch err
                            err isa InterruptException && rethrow()
                            lock(state_lock) do
                                message = sprint(showerror, err)
                                push!(p["errors"], OrderedDict{String,Any}(
                                    "query_index" => j, "error" => message, "at_utc" => utc_timestamp()))
                                first_error[] === nothing && (first_error[] = message)
                                stopped[] = true
                            end
                            break
                        end
                        lock(state_lock) do
                            push!(p["samples"], sample)
                            dirty[] += 1
                            if dirty[] >= checkpoint_every || time()-last_save[] >= checkpoint_seconds
                                s = checkpoint!(data, path)
                                dirty[] = 0
                                last_save[] = time()
                                verbose && println("Saved $(s["saved_responses"])/$(s["planned_queries"]) responses; $(s["completed_prompts"])/$(s["prompt_count"]) prompts complete.")
                            end
                        end
                    end
                  end
                catch
                    # Stop other workers on local failures too (e.g. a failed disk checkpoint).
                    stopped[] = true
                    rethrow()
                end
                push!(tasks, task)
            end
        catch
            # If Ctrl-C interrupts the parent's wait, stop scheduling and drain workers
            # before serializing their shared state. A second interrupt may abort this wait.
            stopped[] = true
            for task in tasks
                try
                    wait(task)
                catch
                    # Preserve the original exception; all workers stop at their next boundary.
                end
            end
            rethrow()
        finally
            # @sync waits for in-flight workers; preserve their responses even if another failed.
            data["run_state"] = first_error[] !== nothing ? "failed" :
                collection_summary(data)["complete"] ? "complete" : "interrupted"
            checkpoint!(data, path)
        end
        first_error[] === nothing || error("Collection stopped; saved responses are in $path. " *
            "Fix the API error and rerun the same command to resume. First error: $(first_error[])")
        verbose && println("Complete. All results saved to $path")
        return data["summary"]
    end
end

end # module
