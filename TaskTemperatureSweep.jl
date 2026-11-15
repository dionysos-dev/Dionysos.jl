using Random
using SHA
using Printf

if !isdefined(@__MODULE__, :request_translation) || !isdefined(@__MODULE__, :build_task_specific_translation_prompt) || !isdefined(@__MODULE__, :TASK_PROVIDER_ADAPTER_VERSION)
    include(joinpath(@__DIR__, "taskSpescific_evaluate.jl"))
end

const PILOT_TEMPERATURES = (0.7,)
const PILOT_REPETITIONS = 10
const PILOT_RESULTS_PATH = joinpath(@__DIR__, "Results", "task-specific-temperature-0.7")
const PILOT_DEV_SPLIT = joinpath(@__DIR__, "Benchmark", "Splits", "train.json")
const PILOT_PROVIDER_SUPPORT_VERSION = 1

pilot_mutable(x::AbstractDict) = OrderedDict{String,Any}(String(k) => pilot_mutable(v) for (k,v) in pairs(x))
pilot_mutable(x::AbstractVector) = Any[pilot_mutable(v) for v in x]
pilot_mutable(x) = x
pilot_read(path) = pilot_mutable(JSON3.read(read(path, String)))
pilot_save(path, obj) = save_results(OrderedDict{String,Any}(obj), path)
pilot_hash(path) = bytes2hex(sha256(read(path)))
pilot_ast_size(record) = get(record, "ast_size", nothing) === nothing ?
    ltl_ast_metrics(parse_ltl(String(record["LTL"]))).ast_size : Int(record["ast_size"])

function pilot_reference_error(record)
    try
        atoms = Set(pilot_atoms(parse_ltl(String(record["LTL"]))))
        unknown = setdiff(atoms, Set(TASK_SPECIFIC_ATOMIC_PROPOSITIONS))
        isempty(unknown) || return "Unavailable reference propositions: $(join(sort!(collect(unknown)), ", "))"
        return nothing
    catch err
        err isa InterruptException && rethrow()
        return sprint(showerror, err)
    end
end

"""Select exactly N prompts, one per template, from training templates only.

The existing splits contain abstract formulas: task-specific template_id, not
task-specific record id, refers to their id. Sample across five-node AST-size
bins, then rotate the four rewritten language fields to balance their sources.
"""
function select_development_prompts(dataset, dev_records; prompt_count=40, seed=20261012)
    rng = MersenneTwister(seed)
    dev_ids = Set(Int(r["id"]) for r in dev_records)
    by_template = Dict{Int,Vector{Any}}()
    for record in dataset
        template = Int(record["template_id"])
        template in dev_ids || continue
        isnothing(pilot_reference_error(record)) || continue
        all(f -> haskey(record, f) && !isempty(strip(String(record[f]))), DEFAULT_INPUT_FIELDS) || continue
        push!(get!(by_template, template, Any[]), record)
    end
    length(by_template) >= prompt_count || throw(ArgumentError("Not enough distinct development templates for $prompt_count prompts"))
    bins = Dict{Int,Vector{Int}}()
    for template in sort!(collect(keys(by_template)))
        size_bin = (pilot_ast_size(first(by_template[template])) - 1) ÷ 5
        push!(get!(bins, size_bin, Int[]), template)
    end
    bin_keys = sort!(collect(keys(bins)))
    foreach(b -> shuffle!(rng, bins[b]), bin_keys)
    chosen = Int[]
    while length(chosen) < prompt_count
        for b in bin_keys
            isempty(bins[b]) && continue
            push!(chosen, pop!(bins[b]))
            length(chosen) == prompt_count && break
        end
    end
    shuffle!(rng, chosen)
    prompts = Any[]
    for (index, template) in enumerate(chosen)
        candidates = sort!(by_template[template]; by=r -> Int(r["id"]))
        record = rand(rng, candidates)
        field = DEFAULT_INPUT_FIELDS[mod1(index, length(DEFAULT_INPUT_FIELDS))]
        push!(prompts, OrderedDict{String,Any}(
            "prompt_index" => index, "record_id" => record["id"],
            "template_id" => template, "input_field" => field,
            "source_text" => String(record[field]), "ground_truth_ltl" => String(record["LTL"]),
            "ast_size" => pilot_ast_size(record), "label" => get(record, "label", nothing),
            "task_specific_ap_mapping" => pilot_mutable(get(record, "task_specific_ap_mapping", Dict{String,Any}())),
            "translation_prompt" => build_task_specific_translation_prompt(String(record[field]);
                ap_mapping=get(record, "task_specific_ap_mapping", nothing)),
        ))
    end
    return prompts
end

"""Spot filter result: 0 means match; 1 means no match; others are errors."""
function pilot_spot_match(formula, filters...)
    executable = Sys.which("ltlfilt")
    isnothing(executable) && throw(ArgumentError("Spot's ltlfilt is required"))
    output, errors = IOBuffer(), IOBuffer()
    cmd = Cmd([executable, "-f", String(formula), "--count", filters...])
    process = run(pipeline(ignorestatus(cmd); stdout=output, stderr=errors))
    process.exitcode in (0, 1) || error("Spot failed: " * strip(String(take!(errors))))
    return process.exitcode == 0
end

function pilot_atoms(node::LTLAtom)
    node.is_constant && return String[]
    name = node.name
    return [startswith(name, '"') ? String(JSON3.read(name)) : name]
end
pilot_atoms(node::LTLUnary) = pilot_atoms(node.child)
pilot_atoms(node::LTLBinary) = vcat(pilot_atoms(node.left), pilot_atoms(node.right))

"""Validity includes syntax, available propositions, and no truncated output.

Sufficiency requires satisfiability and candidate => reference. False/unsatisfiable
formulas cannot earn sufficiency through vacuous implication.
"""
function score_pilot_sample(response, reference)
    candidate = strip(String(response.formula))
    validation_error = nothing
    try
        response.finish_reason == "length" && throw(ArgumentError("Output truncated at max_tokens"))
        ast = parse_ltl(candidate)
        unknown = setdiff(Set(pilot_atoms(ast)), Set(TASK_SPECIFIC_ATOMIC_PROPOSITIONS))
        isempty(unknown) || throw(ArgumentError("Unavailable propositions: $(join(sort!(collect(unknown)), ", "))"))
    catch err
        err isa InterruptException && rethrow()
        validation_error = sprint(showerror, err)
    end
    valid = isnothing(validation_error)
    # Infrastructure/Spot errors propagate; they are never labelled invalid formulas.
    satisfiable = valid && !are_ltl_equivalent(candidate, "0")
    equivalent = satisfiable && are_ltl_equivalent(candidate, reference)
    sufficient = satisfiable && (equivalent || pilot_spot_match(candidate, "--imply=$reference"))
    usage = pilot_mutable(response.usage)
    tokens = isnothing(usage) ? nothing : get(usage, "completion_tokens", get(usage, "output_tokens", nothing))
    return OrderedDict{String,Any}(
        "raw_output" => response.raw_output, "candidate_formula" => candidate,
        "valid" => valid, "validation_error" => validation_error,
        "equivalent" => equivalent, "satisfiable" => satisfiable,
        "semantically_sufficient" => sufficient,
        "usage" => usage, "output_tokens" => tokens,
        "returned_model" => response.returned_model, "finish_reason" => response.finish_reason,
        "generated_at_utc" => string(now(UTC)),
    )
end

function summarize_pilot_prompt(samples)
    clusters = Any[]
    satisfiable_count = 0
    unsatisfiable_count = 0
    sat_cache = Dict{String,Bool}()
    for sample in samples
        get(sample, "valid", false) || continue
        formula = sample["candidate_formula"]
        sat = get(sample,"satisfiable",nothing)
        isnothing(sat) && (sat=get!(() -> !are_ltl_equivalent(formula,"0"),sat_cache,String(formula)))
        if !sat
            unsatisfiable_count += 1
            continue
        end
        satisfiable_count += 1
        cluster_index = findfirst(c -> formula == c["representative"] ||
            are_ltl_equivalent(formula, c["representative"]), clusters)
        if isnothing(cluster_index)
            push!(clusters, OrderedDict{String,Any}(
                "cluster_id" => length(clusters) + 1, "representative" => formula,
                "frequency" => 0, "sample_indices" => Int[],
                "equivalent" => sample["equivalent"],
                "semantically_sufficient" => sample["semantically_sufficient"],
            ))
            cluster_index = length(clusters)
        end
        cluster = clusters[cluster_index]
        cluster["frequency"] += 1
        push!(cluster["sample_indices"], sample["repetition"])
    end
    sort!(clusters; by=c -> (-c["frequency"], first(c["sample_indices"])))
    modal = isempty(clusters) ? nothing : first(clusters)
    tied = isnothing(modal) ? Any[] : filter(c -> c["frequency"] == modal["frequency"], clusters)
    return OrderedDict{String,Any}(
        "samples" => length(samples),
        "valid_count" => count(s -> get(s, "valid", false), samples),
        "satisfiable_count" => satisfiable_count,
        "unsatisfiable_candidate_count" => unsatisfiable_count,
        "clustering_policy" => "valid_and_satisfiable",
        "correct_count" => count(s -> get(s, "equivalent", false), samples),
        "semantic_clusters" => clusters,
        "distinct_semantic_clusters" => length(clusters),
        "modal_cluster_id" => isnothing(modal) ? nothing : modal["cluster_id"],
        "modal_cluster_correct" => !isnothing(modal) && modal["equivalent"],
        "modal_tie" => length(tied) > 1,
        "any_modal_cluster_correct" => any(c -> c["equivalent"], tied),
        "all_modal_clusters_correct" => !isempty(tied) && all(c -> c["equivalent"], tied),
        "correct_formula_appears" => any(s -> get(s, "equivalent", false), samples),
        "sufficient_formula_appears" => any(s -> get(s, "semantically_sufficient", false), samples),
    )
end

pilot_temperature_label(temperature) = isnothing(temperature) ? "default" : @sprintf("%.1f",temperature)
pilot_detail_path(directory, index, temperature) = joinpath(directory, "details",
    @sprintf("prompt_%03d_T_%s.json", index, pilot_temperature_label(temperature)))

function pilot_detail(manifest, prompt, temperature, path)
    if isfile(path)
        result = pilot_read(path)
        result["manifest_sha256"] == bytes2hex(sha256(JSON3.write(manifest))) || error("Checkpoint manifest mismatch at $path")
        result["prompt"] == prompt || error("Checkpoint prompt mismatch at $path")
        result["temperature"] == temperature || error("Checkpoint temperature mismatch at $path")
        repetitions = [s["repetition"] for s in result["samples"]]
        repetitions == collect(1:length(repetitions)) || error("Invalid sample indices at $path")
        length(repetitions) <= PILOT_REPETITIONS || error("Too many samples at $path")
        return result
    end
    return OrderedDict{String,Any}(
        "manifest_sha256" => bytes2hex(sha256(JSON3.write(manifest))),
        "prompt" => prompt, "temperature" => temperature, "samples" => Any[],
    )
end

function prepare_temperature_pilot(; dataset_path=DEFAULT_DATASET_PATH,
    dev_split_path=PILOT_DEV_SPLIT, results_directory=PILOT_RESULTS_PATH,
    prompt_count=40, seed=20261012, model=DEFAULT_MODEL, provider=DEFAULT_PROVIDER,
    temperatures=collect(PILOT_TEMPERATURES), include_temperature=true,
    api_url=provider==:anthropic ? "https://api.anthropic.com/v1/messages" :
        provider==:mistral ? "https://api.mistral.ai/v1/chat/completions" : DEFAULT_API_URL,
    api_key_env=provider==:anthropic ? "CLAUDE_API_KEY" : provider==:mistral ? "MISTRAL_API_KEY" : DEFAULT_API_KEY_ENV,
    thinking_type=provider==:anthropic && startswith(model,"claude-sonnet-5-5") ? "between_tools" :
        provider==:mistral ? "none" : "disabled",
    reasoning_effort=provider==:mistral ? "none" : DEFAULT_REASONING_EFFORT)
    provider in (:deepseek,:anthropic,:mistral) || throw(ArgumentError("Pilot provider must be :deepseek, :anthropic, or :mistral"))
    if provider==:mistral && startswith(model,"mistral-large-4")
        reasoning_effort in ("none","high") || throw(ArgumentError("Mistral Large 4 reasoning_effort must be none or high"))
        thinking_type==reasoning_effort || throw(ArgumentError("Mistral thinking_type must match reasoning_effort"))
    end
    temps = collect(temperatures)
    isempty(temps) && throw(ArgumentError("At least one sampling setting is required"))
    length(unique(temps))==length(temps) || throw(ArgumentError("Temperature settings must be distinct"))
    if include_temperature
        all(t -> t isa Real && isfinite(t) && 0<=t<=(provider==:anthropic ? 1 : 2),temps) ||
            throw(ArgumentError("Invalid temperatures for $provider"))
        foreach(t -> validate_translation_sampling(provider,model,true,Float64(t)),temps)
    else
        temps==[nothing] || throw(ArgumentError("Use temperatures=[nothing] for omitted/provider-default temperature; numeric labels would misrepresent sampling"))
    end
    30 <= prompt_count <= 50 || throw(ArgumentError("Use 30–50 development prompts"))
    isfile(dataset_path) || throw(ArgumentError("Dataset not found: $dataset_path"))
    isfile(dev_split_path) || throw(ArgumentError("Development source split not found: $dev_split_path"))
    manifest_path = joinpath(results_directory, "manifest.json")
    identity = OrderedDict{String,Any}(
        "schema_version" => 2,
        "dataset_path" => abspath(dataset_path), "dataset_sha256" => pilot_hash(dataset_path),
        "development_split_path" => abspath(dev_split_path), "development_split_sha256" => pilot_hash(dev_split_path),
        "prompt_count" => prompt_count, "seed" => seed, "model" => model,
        "temperatures" => temps, "repetitions" => PILOT_REPETITIONS,
        "provider" => string(provider), "include_temperature" => include_temperature,
        "api_key_env" => api_key_env,
        "max_tokens" => DEFAULT_MAX_TOKENS, "thinking_type" => thinking_type,
        "api_url" => api_url,
        "system_prompt" => TASK_TRANSLATION_SYSTEM_PROMPT,
        "prompt_version" => TASK_SPECIFIC_PROMPT_VERSION,
        "atomic_propositions" => copy(TASK_SPECIFIC_ATOMIC_PROPOSITIONS),
        "source_label_resolution" => "Legacy prop_N labels are explicitly mapped to colors in the target prompt using task_specific_ap_mapping",
        "translation_prompt_sha256" => bytes2hex(sha256(build_task_specific_translation_prompt("{requirement}"))),
        "selection" => "Round-robin AST-size bins; one prompt per training template; balanced rewritten fields",
        "sufficiency" => "satisfiable(candidate) AND candidate implies reference; unconstrained infinite traces",
        "modal_tie_break" => "Earliest sampled member among equally frequent valid semantic clusters",
        "validity" => "Parsed LTL over allowed color propositions; finish_reason must not be length; wrappers normalized",
    )
    # Leave existing DeepSeek/Claude manifest identities unchanged on resume.
    provider==:mistral && (identity["reasoning_effort"]=reasoning_effort)
    if isfile(manifest_path)
        manifest = pilot_read(manifest_path)
        legacy_defaults = Dict("provider"=>"deepseek","include_temperature"=>true,"api_key_env"=>DEFAULT_API_KEY_ENV)
        for (key, value) in identity
            get(manifest,key,get(legacy_defaults,key,nothing)) == value || error("Experiment configuration changed ($key); choose a new --results directory")
        end
        return manifest
    end
    dataset = load_dataset(dataset_path)
    prompts = select_development_prompts(dataset, pilot_read(dev_split_path);
        prompt_count=prompt_count, seed=seed)
    identity["excluded_reference_records"] = Any[
        Dict("record_id"=>r["id"], "template_id"=>r["template_id"], "reason"=>pilot_reference_error(r))
        for r in dataset if !isnothing(pilot_reference_error(r))]
    identity["created_at_utc"] = string(now(UTC))
    identity["prompts"] = prompts
    pilot_save(manifest_path, identity)
    return identity
end

function write_pilot_summary(results_directory, manifest)
    rows = Any[]
    prompt_rows = Any[]
    for temperature in manifest["temperatures"]
        samples, metrics = Any[], Any[]
        for prompt in manifest["prompts"]
            path = pilot_detail_path(results_directory, prompt["prompt_index"], temperature)
            isfile(path) || continue
            result = pilot_detail(manifest, prompt, temperature, path)
            append!(samples, result["samples"])
            length(result["samples"]) == PILOT_REPETITIONS || continue
            stats = summarize_pilot_prompt(result["samples"])
            push!(metrics, stats)
            push!(prompt_rows, merge(OrderedDict{String,Any}(
                "prompt_index" => prompt["prompt_index"], "record_id" => prompt["record_id"],
                "input_field" => prompt["input_field"], "temperature" => temperature), stats))
        end
        n, p = length(samples), length(metrics)
        rate(key) = p == 0 ? nothing : count(m -> m[key], metrics) / p
        tokens = [s["output_tokens"] for s in samples if !isnothing(get(s, "output_tokens", nothing))]
        valid_count = count(s -> get(s, "valid", false), samples)
        correct_count = count(s -> get(s, "equivalent", false), samples)
        push!(rows, OrderedDict{String,Any}(
            "temperature" => temperature, "completed_prompts" => p,
            "planned_prompts" => manifest["prompt_count"], "samples" => n,
            "planned_samples" => manifest["prompt_count"] * PILOT_REPETITIONS,
            "valid_ltl_rate" => n == 0 ? nothing : valid_count/n,
            "mean_distinct_semantic_clusters" => p == 0 ? nothing : sum(m["distinct_semantic_clusters"] for m in metrics)/p,
            "per_sample_translation_accuracy" => n == 0 ? nothing : correct_count/n,
            "accuracy_among_valid" => valid_count == 0 ? nothing : correct_count/valid_count,
            "modal_cluster_accuracy" => rate("modal_cluster_correct"),
            "modal_tie_rate" => rate("modal_tie"),
            "modal_accuracy_tie_lower" => rate("all_modal_clusters_correct"),
            "modal_accuracy_tie_upper" => rate("any_modal_cluster_correct"),
            "correct_formula_appears_rate" => rate("correct_formula_appears"),
            "sufficient_formula_appears_rate" => rate("sufficient_formula_appears"),
            "mean_output_tokens" => isempty(tokens) ? nothing : sum(tokens)/length(tokens),
            "total_output_tokens" => sum(tokens; init=0),
            "samples_with_token_usage" => length(tokens),
            "api_error_count" => count(s -> get(s, "status", "ok") == "api_error", samples),
        ))
    end
    summary = OrderedDict{String,Any}(
        "clustering_policy" => "valid_and_satisfiable",
        "updated_at_utc" => string(now(UTC)), "model" => manifest["model"],
        "provider" => get(manifest,"provider","deepseek"),
        "include_temperature" => get(manifest,"include_temperature",true),
        "complete" => all(r -> r["completed_prompts"] == manifest["prompt_count"], rows),
        "denominators" => "Sample rates include invalid outputs and exhausted API errors; prompt rates use complete ten-sample groups; clusters exclude invalid and individually unsatisfiable outputs; token mean uses reported usage only",
        "temperature_results" => rows, "per_prompt_results" => prompt_rows,
    )
    pilot_save(joinpath(results_directory, "summary.json"), summary)
    fields = collect(keys(first(rows)))
    csv_path = joinpath(results_directory, "summary.csv")
    open(csv_path, "w") do io
        println(io, join(fields, ','))
        for row in rows
            println(io, join([f=="temperature" ? pilot_temperature_label(row[f]) :
                (isnothing(row[f]) ? "" : string(row[f])) for f in fields], ','))
        end
    end
    report_path = joinpath(results_directory, "comparison.md")
    open(report_path, "w") do io
        println(io, "Temperature pilot — $(manifest["model"])")
        println(io, "\n$(manifest["prompt_count"]) fixed development prompts × $(manifest["repetitions"]) samples × $(length(manifest["temperatures"])) temperature(s). Complete: $(summary["complete"]).\n")
        if !get(manifest,"include_temperature",true)
            println(io,"Temperature is provider-default: the temperature parameter is omitted from requests. This is a single sampling setting, not a temperature-tuning sweep. Thinking setting: $(manifest["thinking_type"]).\n")
        end
        if haskey(manifest,"reasoning_effort")
            println(io,"Reasoning effort is fixed at $(manifest["reasoning_effort"]) across all temperatures.\n")
        end
        println(io, "| T | Complete prompts | Valid LTL | Mean clusters | Sample accuracy | Modal accuracy | Correct appears | Sufficient appears | Mean output tokens |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in rows
            pct(k) = isnothing(row[k]) ? "—" : @sprintf("%.1f%%", row[k]*100)
            num(k) = isnothing(row[k]) ? "—" : @sprintf("%.2f", row[k])
            println(io, "| $(pilot_temperature_label(row["temperature"])) | $(row["completed_prompts"]) | $(pct("valid_ltl_rate")) | $(num("mean_distinct_semantic_clusters")) | $(pct("per_sample_translation_accuracy")) | $(pct("modal_cluster_accuracy")) | $(pct("correct_formula_appears_rate")) | $(pct("sufficient_formula_appears_rate")) | $(num("mean_output_tokens")) |")
        end
        println(io, "\nSufficiency means satisfiable candidate ⇒ reference under unconstrained infinite-trace LTL. Ties use the cluster containing the earliest sample; JSON/CSV also report tie rate and lower/upper modal accuracy. All-invalid groups count as incorrect and have zero clusters.")
        println(io, "\nCandidates are filtered for individual satisfiability before semantic clustering and modal selection. Output-token cost still includes excluded candidates; per-prompt JSON reports their counts.")
        println(io, "\nToken cost is measured in API-reported output tokens, including invalid generations. Missing usage is excluded from that average and its coverage is reported. Sample accuracy includes every scheduled sample. API failures are recorded separately. The pilot is for development selection; it does not estimate held-out test performance.")
        println(io, "\nExcluded $(length(manifest["excluded_reference_records"])) benchmark records with malformed references or propositions outside the five allowed colors. Their IDs and reasons are saved in manifest.json; the source benchmark is unchanged.")
    end
    return summary
end

function run_temperature_pilot(; workers=6, request_fn=request_translation, kwargs...)
    workers >= 1 || throw(ArgumentError("workers must be positive"))
    manifest = prepare_temperature_pilot(; kwargs...)
    directory = get(kwargs, :results_directory, PILOT_RESULTS_PATH)
    Sys.which("ltlfilt") === nothing && error("Install Spot before running")
    # Check every selected reference before incurring any API cost.
    for prompt in manifest["prompts"]
        pilot_spot_match(prompt["ground_truth_ltl"])
        actual_prompt = build_task_specific_translation_prompt(prompt["source_text"];
            ap_mapping=prompt["task_specific_ap_mapping"])
        actual_prompt == prompt["translation_prompt"] || error("Saved translation prompt changed; use a new results directory")
    end
    request_fn === request_translation && get_api_key(get(manifest,"api_key_env",DEFAULT_API_KEY_ENV))
    jobs = [(p,t) for p in manifest["prompts"] for t in manifest["temperatures"]]
    # Seeded interleaving prevents all of one temperature being run first.
    shuffle!(MersenneTwister(manifest["seed"] + 1), jobs)
    queue = Channel{Any}(length(jobs))
    foreach(job -> put!(queue, job), jobs)
    close(queue)
    completed = Ref(0)
    try
        @sync for _ in 1:min(workers, length(jobs))
            @async for (prompt, temperature) in queue
                path = pilot_detail_path(directory, prompt["prompt_index"], temperature)
                result = pilot_detail(manifest, prompt, temperature, path)
                for repetition in (length(result["samples"])+1):PILOT_REPETITIONS
                    pending_path = path * ".pending.json"
                    response = try
                        pending = isfile(pending_path) ? pilot_read(pending_path) : nothing
                        if !isnothing(pending) && pending["repetition"] == repetition
                            saved = pending["response"]
                            (formula=saved["formula"], raw_output=saved["raw_output"],
                                usage=saved["usage"], returned_model=saved["returned_model"],
                                finish_reason=saved["finish_reason"])
                        else
                            request_fn(prompt["source_text"]; provider=Symbol(get(manifest,"provider","deepseek")),
                                model=manifest["model"], temperature=isnothing(temperature) ? 1.0 : Float64(temperature),
                                include_temperature=get(manifest,"include_temperature",true),
                                api_url=manifest["api_url"],api_key_env=get(manifest,"api_key_env",DEFAULT_API_KEY_ENV),
                                max_tokens=manifest["max_tokens"], deepseek_thinking_type=manifest["thinking_type"],
                                reasoning_effort=get(manifest,"reasoning_effort",DEFAULT_REASONING_EFFORT),
                                anthropic_thinking_type=manifest["thinking_type"],
                                ap_mapping=prompt["task_specific_ap_mapping"],
                                return_metadata=true)
                        end
                    catch err
                        err isa InterruptException && rethrow()
                        err isa Downloads.RequestError && rethrow()
                        # Stop the run on authentication, model, or balance errors.
                        message = sprint(showerror, err)
                        occursin(r"status (400|401|402|403|404)", message) && rethrow()
                        sample = OrderedDict{String,Any}("repetition" => repetition,
                            "status" => "api_error", "error" => message,
                            "valid" => false, "equivalent" => false,
                            "semantically_sufficient" => false, "output_tokens" => nothing,
                            "generated_at_utc" => string(now(UTC)))
                        push!(result["samples"], sample)
                        pilot_save(path, result)
                        continue
                    end
                    # Persist the response before semantic analysis, so Spot failures
                    # do not require paying for the same sample again on resume.
                    pending = OrderedDict{String,Any}("repetition" => repetition,
                        "response" => Dict(string(k)=>pilot_mutable(v) for (k,v) in pairs(response)))
                    pilot_save(pending_path, pending)
                    sample = score_pilot_sample(response, prompt["ground_truth_ltl"])
                    sample["repetition"] = repetition
                    sample["status"] = "ok"
                    push!(result["samples"], sample)
                    pilot_save(path, result)
                end
                completed[] += 1
                println("Completed $(completed[])/$(length(jobs)): prompt=$(prompt["prompt_index"]) T=$(pilot_temperature_label(temperature))")
                # One process owns this directory; Julia async tasks share one thread.
                completed[] % 6 == 0 && write_pilot_summary(directory, manifest)
            end
        end
    finally
        write_pilot_summary(directory, manifest)
    end
    return write_pilot_summary(directory, manifest)
end

function temperature_pilot_main(args=ARGS)
    options = Dict{Symbol,Any}()
    action, workers = "--run", 6
    names = Dict("--prompts" => (:prompt_count, Int), "--seed" => (:seed, Int),
        "--model" => (:model, String), "--results" => (:results_directory, String),
        "--dataset" => (:dataset_path, String), "--dev-split" => (:dev_split_path, String))
    index = 1
    while index <= length(args)
        arg = args[index]
        if arg == "--help"
            println("julia --project=. TaskTemperatureSweep.jl [--prepare|--run|--summarize] [--prompts 40] [--seed 20261012] [--workers 6] [--model deepseek-v4-pro] [--results DIR] [--dataset FILE] [--dev-split FILE]")
            return
        elseif arg in ("--prepare", "--run", "--summarize")
            action = arg
        elseif arg == "--workers" || haskey(names, arg)
            index += 1
            index <= length(args) || error("Missing value for $arg")
            if arg == "--workers"
                workers = parse(Int, args[index])
            else
                key, typ = names[arg]
                options[key] = typ == String ? args[index] : parse(typ, args[index])
            end
        else
            error("Unknown argument: $arg")
        end
        index += 1
    end
    directory = get(options, :results_directory, PILOT_RESULTS_PATH)
    if action == "--run"
        run_temperature_pilot(; workers=workers, options...)
    else
        manifest = prepare_temperature_pilot(; options...)
        write_pilot_summary(directory, manifest)
        println("Prepared $(manifest["prompt_count"]) prompts; $(manifest["prompt_count"]*manifest["repetitions"]*length(manifest["temperatures"])) planned outputs; sampling settings: $(join(pilot_temperature_label.(manifest["temperatures"]),", ")).")
    end
    println("Comparison: ", joinpath(directory, "comparison.md"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    temperature_pilot_main()
end
