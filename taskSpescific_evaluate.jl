

if !isdefined(@__MODULE__, :are_ltl_equivalent)
    include(joinpath(@__DIR__, "LTL.jl"))
end
are_equivalent(a, b) = are_ltl_equivalent(a, b)

using JSON3
using OrderedCollections
using Downloads
using Dates

# =================================================================================================
# CONFIGURATION
# Change only these values to evaluate a different model and save to a different JSON file.
# =================================================================================================

const DEFAULT_DATASET_PATH = joinpath(@__DIR__, "Benchmark", "Task_Spescific_Dataset.json")


const DEFAULT_RESULTS_PATH = joinpath(
    @__DIR__,
    "Results",
    "DeepSeekV4Pro_task_specific_rewritten.json",
)
const DEFAULT_PROMPT_SETTING = "fewshot"

const DEFAULT_PROVIDER = :deepseek
const DEFAULT_MODEL = "deepseek-v4-pro"
const DEFAULT_API_URL = "https://api.deepseek.com/chat/completions"
const DEFAULT_API_KEY_ENV = "DEEPSEEK_API_KEY"
const DEFAULT_INCLUDE_TEMPERATURE = true
const DEFAULT_TEMPERATURE = 0.7
const DEFAULT_ANTHROPIC_VERSION = "2023-06-01"
const DEFAULT_MAX_TOKENS = 1024
const DEFAULT_REASONING_EFFORT = "disabled"
const DEFAULT_DEEPSEEK_THINKING_TYPE = "disabled"

const DEFAULT_INPUT_FIELDS = [
    "natural_paraphrase_gpt54mini_rewritten",
    "paraphrase_gpt5.4-mini_gpt54mini_rewritten",
    "paraphrase_gemini-2.5-flash_gpt54mini_rewritten",
    "paraphrase_deepseek_gpt54mini_rewritten",
]

const TASK_SPECIFIC_ATOMIC_PROPOSITIONS = [
    "blue",
    "brown",
    "green",
    "purple",
    "yellow",
]
const TASK_SPECIFIC_PROMPT_VERSION = 2
const TASK_TRANSLATION_SYSTEM_PROMPT = "You translate natural language statements into LTL while preserving exact meaning."
const TASK_PROVIDER_ADAPTER_VERSION = 1

function validate_translation_sampling(provider, model, include_temperature, temperature)
    if provider==:anthropic && startswith(model,"claude-sonnet-5-5") && include_temperature && temperature!=1.0
        throw(ArgumentError("Claude Sonnet 5.5 rejects non-default temperatures. Use include_temperature=false and provider-default sampling; do not label omitted temperatures as a tuning sweep."))
    end
end

# =================================================================================================
# I/O helpers
# =================================================================================================

function load_dataset(dataset_path::String = DEFAULT_DATASET_PATH)
    if !isfile(dataset_path)
        throw(ArgumentError("Dataset file not found: $(dataset_path)"))
    end
    data = JSON3.read(read(dataset_path, String))
    return [OrderedDict{String,Any}(String(k) => v for (k, v) in pairs(record)) for record in data]
end

function ensure_parent_directory(path::String)
    parent = dirname(path)
    isdir(parent) || mkpath(parent)
end

function save_results(results_obj::OrderedDict{String,Any}, results_path::String)
    ensure_parent_directory(results_path)
    parent = dirname(results_path)
    temp_path, temp_io = mktemp(parent)
    try
        JSON3.pretty(temp_io, results_obj)
        write(temp_io, "\n")
        close(temp_io)
        mv(temp_path, results_path; force=true)
    catch
        isopen(temp_io) && close(temp_io)
        isfile(temp_path) && rm(temp_path)
        rethrow()
    end
end

function initial_results_object(
    dataset_path::String;
    provider::Symbol,
    model::String,
    input_fields::Vector{String},
)
    return OrderedDict(
        "dataset_path" => dataset_path,
        "translation_provider" => string(provider),
        "translation_model" => model,
        "prompt_setting" => DEFAULT_PROMPT_SETTING,
        "input_fields" => input_fields,
        "created_at" => string(now()),
        "updated_at" => string(now()),
        "results" => OrderedDict[],
    )
end

function load_results(
    results_path::String;
    dataset_path::String,
    provider::Symbol,
    model::String,
    input_fields::Vector{String},
)
    if !isfile(results_path)
        return initial_results_object(dataset_path; provider=provider, model=model, input_fields=input_fields)
    end

    content = strip(read(results_path, String))
    isempty(content) && return initial_results_object(dataset_path; provider=provider, model=model, input_fields=input_fields)

    parsed = JSON3.read(content)
    results_obj = OrderedDict{String,Any}(String(k) => v for (k, v) in pairs(parsed))

    if haskey(results_obj, "results")
        materialized_results = OrderedDict{String,Any}[]
        for entry in results_obj["results"]
            push!(materialized_results, OrderedDict{String,Any}(String(k) => v for (k, v) in pairs(entry)))
        end
        results_obj["results"] = materialized_results
    else
        results_obj["results"] = OrderedDict[]
    end

    return results_obj
end

# =================================================================================================
# API helpers
# =================================================================================================

function get_api_key(env_name::String)
    env_path = joinpath(@__DIR__, ".env")
    if isfile(env_path)
        for raw_line in eachline(env_path)
            line = strip(raw_line)
            (isempty(line) || startswith(line, '#') || !occursin('=', line)) && continue
            startswith(line, "export ") && (line = strip(line[8:end]))
            key, value = strip.(split(line, '='; limit=2))
            if length(value) >= 2 && first(value) == last(value) && first(value) in ('\"', '\'')
                value = value[2:end-1]
            end
            isempty(key) || haskey(ENV, key) || (ENV[key] = value)
        end
    end

    api_key = get(ENV, env_name, "")
    isempty(api_key) && throw(ArgumentError(
        "$(env_name) is not set. Put it in a local .env file or export it in your shell before running this script."
    ))
    return api_key
end

"""POST JSON using Julia's standard library; retry only transient failures."""
function post_json(api_url, headers, body; max_retries=4, timeout=120)
    downloader = Downloads.Downloader()
    downloader.easy_hook = (easy, _) -> Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_LOW_SPEED_TIME, 0)
    for attempt in 1:max_retries
        output = IOBuffer()
        response = try
            Downloads.request(api_url; method="POST", headers=headers,
                input=IOBuffer(body), output=output, throw=false, timeout=timeout,
                downloader=downloader)
        catch err
            err isa InterruptException && rethrow()
            attempt == max_retries && rethrow()
            sleep(min(8, 2^(attempt - 1)))
            continue
        end
        if response isa Downloads.RequestError
            attempt == max_retries && throw(response)
            sleep(min(8, 2^(attempt - 1)))
            continue
        end
        result = (status=response.status, body=take!(output))
        if (response.status == 429 || response.status >= 500) && attempt < max_retries
            sleep(min(8, 2^(attempt - 1)))
            continue
        end
        return result
    end
    error("API retry loop exhausted")
end

function translation_response(parsed, raw_output; return_metadata=false)
    formula = normalize_formula_text(raw_output)
    return_metadata || return formula
    choice = haskey(parsed, "choices") && !isempty(parsed["choices"]) ? parsed["choices"][1] : nothing
    finish_reason = isnothing(choice) ? get(parsed,"stop_reason",nothing) : get(choice,"finish_reason",nothing)
    finish_reason=="max_tokens" && (finish_reason="length")
    return (
        formula=formula,
        raw_output=String(raw_output),
        usage=get(parsed, "usage", nothing),
        returned_model=get(parsed, "model", nothing),
        finish_reason=finish_reason,
    )
end

"""Build the task-specific prompt, including AP meanings and source-label aliases.

Every few-shot formula uses colors from the target vocabulary. Legacy prop_N
labels occurring in the requirement are resolved using the record's saved
mapping. No reference formula is exposed to the model.
"""
function build_task_specific_translation_prompt(sentence::AbstractString; ap_mapping=nothing)
    available_atomic_propositions = "AP = {" * join(TASK_SPECIFIC_ATOMIC_PROPOSITIONS, ", ") * "}"
    legacy_labels = sort!(unique([m.match for m in eachmatch(r"prop_\d+(?![A-Za-z0-9_])", sentence)]))
    aliases = String[]
    for label in legacy_labels
        isnothing(ap_mapping) && throw(ArgumentError("Requirement contains $label; provide its task_specific_ap_mapping"))
        haskey(ap_mapping, label) || throw(ArgumentError("Missing task_specific_ap_mapping entry for $label"))
        color = ap_mapping[label]
        color in TASK_SPECIFIC_ATOMIC_PROPOSITIONS || throw(ArgumentError("$label maps to an unavailable task proposition: $color"))
        push!(aliases, "- `$(label)` means `$(color)` (the robot is in the $(color) region).")
    end
    alias_section = isempty(aliases) ? String[] : vcat([
        "Legacy labels appearing in this requirement:",
        "These labels are aliases for the color propositions; output their color names.",
    ], aliases, [""])

    return join(vcat([
        "You translate one natural-language requirement into a Linear Temporal Logic (LTL) formula.",
        "",
        "Task atomic propositions:",
        available_atomic_propositions,
        "Each color proposition is true exactly when the robot is in that color's region at the current time.",
        "Use only these five color names as atomic propositions in the output.",
        "",
        "Semantics:",
        "",
        "- Use standard LTL over infinite traces, evaluated from the current time.",
        "- `U` is strong Until: `(a U b)` requires `b` to eventually hold, while `a` holds at every preceding time.",
        "- “Whenever,” “every time,” and “globally” introduce a recurring global requirement.",
        "- “From then on” means that the consequent must hold globally from the triggering time.",
        "- Do not add a global operator unless the requirement expresses a global or recurring condition.",
        "",
        "Allowed syntax:",
        "",
        "- Atomic propositions: `blue`, `brown`, `green`, `purple`, `yellow`.",
        "- Negation: `!`",
        "- Conjunction: `&`",
        "- Disjunction: `|`",
        "- Implication: `->`",
        "- Equivalence: `<->`",
        "- Next: `X`",
        "- Eventually: `F`",
        "- Globally: `G`",
        "- Until: `U`",
        "",
        "Translation rules:",
        "",
        "- Preserve the complete logical and temporal meaning.",
        "- Use the exact lowercase color names listed in AP.",
        "- Do not introduce new propositions.",
        "- Do not introduce domain assumptions, such as mutual exclusion, reachability, or ordering, unless explicitly stated.",
        "- Do not use every available proposition unless the requirement needs it.",
        "- Use parentheses to make operator scope unambiguous.",
        "- Before answering, silently check temporal scope, ordering, eventuality, and whether each Until obligation is strong.",
        "",
        "Output rules:",
        "",
        "- Output exactly one LTL formula.",
        "- Do not output explanations, labels, comments, or code fences.",
        "",
        "Examples:",
        "",
        "Natural-language requirement:",
        "Whenever the robot is in the blue region, it must remain in the yellow region from then on.",
        "",
        "LTL:",
        "G((blue) -> G(yellow))",
        "",
        "Natural-language requirement:",
        "Whenever the robot is in the blue region, it must stay in the green region until it reaches the yellow region.",
        "",
        "LTL:",
        "G((blue) -> ((green) U (yellow)))",
        "",
        "Natural-language requirement:",
        "Eventually reach the blue region, and from that time eventually reach the yellow region.",
        "",
        "LTL:",
        "F((blue) & F(yellow))",
        "",
        "Natural-language requirement:",
        "The robot must never enter the purple region.",
        "",
        "LTL:",
        "G(!(purple))",
        "",
        "Natural-language requirement:",
        "Stay in the blue region until reaching the brown region.",
        "",
        "LTL:",
        "(blue) U (brown)",
        "",
        "Natural-language requirement:",
        "At the next step, the robot must be in the green region.",
        "",
        "LTL:",
        "X(green)",
        "",
        "Target input:",
        "",
        "Available atomic proposition set:",
        available_atomic_propositions,
        "",
    ], alias_section, [
        "Natural-language requirement:",
        "<requirement>",
        String(sentence),
        "</requirement>",
        "",
        "LTL:",
    ]), "\n")
end

# Keep the existing public helper, while API requests use the uniquely named
# builder so KRecall.jl's generic build_translation_prompt cannot replace it.
build_translation_prompt(sentence::AbstractString; kwargs...) =
    build_task_specific_translation_prompt(sentence; kwargs...)

function extract_openai_output_text(response_json)
    if haskey(response_json, :output_text)
        return String(response_json[:output_text])
    elseif haskey(response_json, "output_text")
        return String(response_json["output_text"])
    end

    key = haskey(response_json, :output) ? :output : (haskey(response_json, "output") ? "output" : nothing)
    isnothing(key) && throw(ErrorException("OpenAI response did not contain `output_text` or `output`."))

    outputs = response_json[key]
    for item in outputs
        if (haskey(item, :type) && item[:type] == "message") || (haskey(item, "type") && item["type"] == "message")
            content_key = haskey(item, :content) ? :content : "content"
            for part in item[content_key]
                part_type = haskey(part, :type) ? part[:type] : part["type"]
                if part_type in ("output_text", "text")
                    text_key = haskey(part, :text) ? :text : "text"
                    return String(part[text_key])
                end
            end
        end
    end

    throw(ErrorException("Could not extract text from OpenAI response."))
end

function extract_gemini_output_text(response_json)
    candidates = if haskey(response_json, :candidates)
        response_json[:candidates]
    elseif haskey(response_json, "candidates")
        response_json["candidates"]
    else
        throw(ErrorException("Gemini response did not contain `candidates`."))
    end

    isempty(candidates) && throw(ErrorException("Gemini response contained no candidates."))
    first_candidate = candidates[1]
    content = haskey(first_candidate, :content) ? first_candidate[:content] : first_candidate["content"]
    parts = haskey(content, :parts) ? content[:parts] : content["parts"]

    text_chunks = String[]
    for part in parts
        if haskey(part, :text)
            push!(text_chunks, String(part[:text]))
        elseif haskey(part, "text")
            push!(text_chunks, String(part["text"]))
        end
    end

    isempty(text_chunks) && throw(ErrorException("Gemini response contained no text parts."))
    return strip(join(text_chunks, "\n"))
end

function extract_deepseek_output_text(response_json)
    choices = if haskey(response_json, :choices)
        response_json[:choices]
    elseif haskey(response_json, "choices")
        response_json["choices"]
    else
        throw(ErrorException("DeepSeek response did not contain `choices`."))
    end

    isempty(choices) && throw(ErrorException("DeepSeek response did not contain no choices."))
    first_choice = choices[1]
    message = haskey(first_choice, :message) ? first_choice[:message] : first_choice["message"]

    if haskey(message, :content)
        return strip(String(message[:content]))
    elseif haskey(message, "content")
        return strip(String(message["content"]))
    else
        throw(ErrorException("DeepSeek response did not contain `message.content`."))
    end
end

function extract_anthropic_output_text(response_json)
    content = if haskey(response_json, :content)
        response_json[:content]
    elseif haskey(response_json, "content")
        response_json["content"]
    else
        throw(ErrorException("Anthropic response did not contain `content`."))
    end

    text_chunks = String[]
    for part in content
        part_type = haskey(part, :type) ? String(part[:type]) : String(part["type"])
        if part_type == "text"
            if haskey(part, :text)
                push!(text_chunks, String(part[:text]))
            elseif haskey(part, "text")
                push!(text_chunks, String(part["text"]))
            end
        end
    end

    isempty(text_chunks) && throw(ErrorException("Anthropic response contained no text content."))
    return strip(join(text_chunks, "\n"))
end

function extract_mistral_output_text(response_json)
    choices = if haskey(response_json, :choices)
        response_json[:choices]
    elseif haskey(response_json, "choices")
        response_json["choices"]
    else
        throw(ErrorException("Mistral response did not contain `choices`."))
    end

    isempty(choices) && throw(ErrorException("Mistral response contained no choices."))
    first_choice = choices[1]
    message = haskey(first_choice, :message) ? first_choice[:message] : first_choice["message"]
    content = haskey(message, :content) ? message[:content] : message["content"]

    if content isa AbstractString
        return strip(String(content))
    end

    text_chunks = String[]
    for part in content
        part_type = haskey(part, :type) ? String(part[:type]) : get(part,"type","text")
        part_type=="text" || continue
        if haskey(part, :text)
            push!(text_chunks, String(part[:text]))
        elseif haskey(part, "text")
            push!(text_chunks, String(part["text"]))
        elseif haskey(part, :content)
            push!(text_chunks, String(part[:content]))
        elseif haskey(part, "content")
            push!(text_chunks, String(part["content"]))
        end
    end

    isempty(text_chunks) && throw(ErrorException("Mistral response contained no text content."))
    return strip(join(text_chunks, "\n"))
end

function normalize_formula_text(text::AbstractString)
    s = strip(String(text))

    if startswith(s, "```")
        parts = split(s, "```")
        for part in parts
            cleaned = strip(replace(part, "ltl" => "", "LTL" => ""))
            if !isempty(cleaned)
                s = cleaned
                break
            end
        end
    end

    if startswith(s, "LTL:")
        s = strip(s[5:end])
    elseif startswith(s, "ltl:")
        s = strip(s[5:end])
    end

    return strip(s)
end

function request_translation(
    sentence::AbstractString;
    provider::Symbol = DEFAULT_PROVIDER,
    model::String = DEFAULT_MODEL,
    api_url::String = DEFAULT_API_URL,
    api_key_env::String = DEFAULT_API_KEY_ENV,
    include_temperature::Bool = DEFAULT_INCLUDE_TEMPERATURE,
    temperature::Float64 = DEFAULT_TEMPERATURE,
    anthropic_version::String = DEFAULT_ANTHROPIC_VERSION,
    max_tokens::Int = DEFAULT_MAX_TOKENS,
    reasoning_effort::String = DEFAULT_REASONING_EFFORT,
    deepseek_thinking_type::String = DEFAULT_DEEPSEEK_THINKING_TYPE,
    return_metadata::Bool = false,
    ap_mapping = nothing,
    anthropic_thinking_type::String = startswith(model,"claude-sonnet-5-5") ? "between_tools" : "disabled",
)
    validate_translation_sampling(provider,model,include_temperature,temperature)
    prompt = build_task_specific_translation_prompt(sentence; ap_mapping=ap_mapping)
    api_key = get_api_key(api_key_env)

    if provider == :openai
        body = Dict(
            "model" => model,
            "input" => prompt,
        )

        if include_temperature
            body["temperature"] = temperature
        end
        headers = [
            "Authorization" => "Bearer $(api_key)",
            "Content-Type" => "application/json",
        ]
        response = post_json(api_url, headers, JSON3.write(body))
        response.status == 200 || throw(ErrorException(
            "OpenAI API request failed with status $(response.status): $(String(response.body))"
        ))
        parsed = JSON3.read(String(response.body))
        return translation_response(parsed, extract_openai_output_text(parsed); return_metadata=return_metadata)

    elseif provider == :gemini
        endpoint = endswith(api_url, ":generateContent") ? api_url : "$(api_url)/$(model):generateContent"
        body = Dict(
            "contents" => [
                Dict("parts" => [Dict("text" => prompt)])
            ],
            "generationConfig" => Dict("temperature" => temperature),
        )
        headers = [
            "x-goog-api-key" => api_key,
            "Content-Type" => "application/json",
        ]
        response = post_json(endpoint, headers, JSON3.write(body))
        response.status == 200 || throw(ErrorException(
            "Gemini API request failed with status $(response.status): $(String(response.body))"
        ))
        parsed = JSON3.read(String(response.body))
        return translation_response(parsed, extract_gemini_output_text(parsed); return_metadata=return_metadata)

    elseif provider == :deepseek
        body = Dict(
            "model" => model,
            "messages" => [
                Dict("role" => "system", "content" => TASK_TRANSLATION_SYSTEM_PROMPT),
                Dict("role" => "user", "content" => prompt),
            ],
            "max_tokens" => max_tokens,
            "stream" => false,
            "thinking" => Dict("type" => deepseek_thinking_type),
        )
        include_temperature && (body["temperature"]=temperature)
        headers = [
            "Authorization" => "Bearer $(api_key)",
            "Content-Type" => "application/json",
        ]
        response = post_json(api_url, headers, JSON3.write(body))
        response.status == 200 || throw(ErrorException(
            "DeepSeek API request failed with status $(response.status): $(String(response.body))"
        ))
        parsed = JSON3.read(String(response.body))
        return translation_response(parsed, extract_deepseek_output_text(parsed); return_metadata=return_metadata)

    elseif provider == :mistral
        body = Dict(
            "model" => model,
            "messages" => [
                Dict("role" => "system", "content" => TASK_TRANSLATION_SYSTEM_PROMPT),
                Dict("role" => "user", "content" => prompt),
            ],
            "max_tokens" => max_tokens,
        )
        if include_temperature
            body["temperature"] = temperature
        end
        reasoning_effort != "disabled" && (body["reasoning_effort"]=reasoning_effort)
        headers = [
            "Authorization" => "Bearer $(api_key)",
            "Content-Type" => "application/json",
            "Accept" => "application/json",
        ]
        response = post_json(api_url, headers, JSON3.write(body))
        response.status == 200 || throw(ErrorException(
            "Mistral API request failed with status $(response.status): $(String(response.body))"
        ))
        parsed = JSON3.read(String(response.body))
        return translation_response(parsed, extract_mistral_output_text(parsed); return_metadata=return_metadata)

    elseif provider == :anthropic
        startswith(model,"claude-sonnet-5-5") && anthropic_thinking_type=="disabled" &&
            throw(ArgumentError("Sonnet 5.5 uses thinking.type=between_tools to disable up-front thinking"))
        body = Dict(
            "model" => model,
            "max_tokens" => max_tokens,
            "system" => TASK_TRANSLATION_SYSTEM_PROMPT,
            "thinking" => Dict("type"=>anthropic_thinking_type),
            "messages" => [
                Dict(
                    "role" => "user",
                    "content" => prompt,
                )
            ],
        )
        include_temperature && (body["temperature"]=temperature)
        headers = [
            "x-api-key" => api_key,
            "anthropic-version" => anthropic_version,
            "content-type" => "application/json",
        ]
        response = post_json(api_url, headers, JSON3.write(body))
        response.status == 200 || throw(ErrorException(
            "Anthropic API request failed with status $(response.status): $(String(response.body))"
        ))
        parsed = JSON3.read(String(response.body))
        return translation_response(parsed, extract_anthropic_output_text(parsed); return_metadata=return_metadata)
    else
        throw(ArgumentError("Unsupported provider: $(provider). Use :openai, :gemini, :deepseek, :anthropic, or :mistral."))
    end
end

# =================================================================================================
# Result construction
# =================================================================================================

function paraphrase_source_model(field::AbstractString)
    if field == "natural_paraphrase_gpt54mini_rewritten"
        return "natural_paraphrase_gpt54mini_rewritten"
    elseif field == "paraphrase_gpt5.4-mini_gpt54mini_rewritten"
        return "paraphrase_gpt5.4-mini_gpt54mini_rewritten"
    elseif field == "paraphrase_gemini-2.5-flash_gpt54mini_rewritten"
        return "paraphrase_gemini-2.5-flash_gpt54mini_rewritten"
    elseif field == "paraphrase_deepseek_gpt54mini_rewritten"
        return "paraphrase_deepseek_gpt54mini_rewritten"
    else
        return "unknown"
    end
end

function result_key(record_id, field::AbstractString)
    return "$(record_id)||$(field)"
end

function latest_results_by_key(results_obj::OrderedDict{String,Any})
    latest = Dict{String,OrderedDict{String,Any}}()
    for entry_any in results_obj["results"]
        entry = OrderedDict{String,Any}(String(k) => v for (k, v) in pairs(entry_any))
        if haskey(entry, "record_id") && haskey(entry, "input_field")
            rid = entry["record_id"]
            field = String(entry["input_field"])
            latest[result_key(rid, field)] = entry
        end
    end
    return latest
end

function successful_result_keys(results_obj::OrderedDict{String,Any})
    keys_set = Set{String}()
    for (key, entry) in latest_results_by_key(results_obj)
        if haskey(entry, "status") && String(entry["status"]) == "ok"
            push!(keys_set, key)
        end
    end
    return keys_set
end

function append_result!(results_obj::OrderedDict{String,Any}, entry::OrderedDict{String,Any})
    new_key = result_key(entry["record_id"], String(entry["input_field"]))
    retained = OrderedDict{String,Any}[]

    for existing_any in results_obj["results"]
        existing = OrderedDict{String,Any}(
            String(k) => v for (k, v) in pairs(existing_any)
        )
        if haskey(existing, "record_id") && haskey(existing, "input_field")
            existing_key = result_key(
                existing["record_id"],
                String(existing["input_field"]),
            )
            existing_key == new_key && continue
        end
        push!(retained, existing)
    end

    push!(retained, entry)
    results_obj["results"] = retained
    results_obj["updated_at"] = string(now())
    return results_obj
end

function make_result_entry(
    record::OrderedDict{String,Any},
    input_field::String,
    source_text::AbstractString,
    predicted_ltl::Union{Nothing,AbstractString},
    equivalent::Union{Nothing,Bool},
    status::String;
    translation_provider::Symbol,
    translation_model::String,
    error_message::Union{Nothing,String} = nothing,
)
    entry = OrderedDict{String,Any}(
        "record_id" => get(record, "id", nothing),
        "template_id" => get(record, "template_id", nothing),
        "label" => get(record, "label", nothing),
        "input_field" => input_field,
        "source_model" => paraphrase_source_model(input_field),
        "translation_provider" => string(translation_provider),
        "translation_model" => translation_model,
        "prompt_setting" => DEFAULT_PROMPT_SETTING,
        "source_text" => String(source_text),
        "original_ltl" => String(get(record, "LTL", "")),
        "predicted_ltl" => isnothing(predicted_ltl) ? nothing : String(predicted_ltl),
        "equivalent" => equivalent,
        "status" => status,
        "timestamp" => string(now()),
    )

    !isnothing(error_message) && (entry["error"] = error_message)
    return entry
end

# =================================================================================================
# Evaluation loop
# =================================================================================================

function evaluate_model_on_benchmark(
    dataset_path::String = DEFAULT_DATASET_PATH;
    results_path::String = DEFAULT_RESULTS_PATH,
    input_fields::Vector{String} = DEFAULT_INPUT_FIELDS,
    provider::Symbol = DEFAULT_PROVIDER,
    model::String = DEFAULT_MODEL,
    api_url::String = DEFAULT_API_URL,
    api_key_env::String = DEFAULT_API_KEY_ENV,
    include_temperature::Bool = DEFAULT_INCLUDE_TEMPERATURE,
    temperature::Float64 = DEFAULT_TEMPERATURE,
    anthropic_version::String = DEFAULT_ANTHROPIC_VERSION,
    max_tokens::Int = DEFAULT_MAX_TOKENS,
    reasoning_effort::String = DEFAULT_REASONING_EFFORT,
    deepseek_thinking_type::String = DEFAULT_DEEPSEEK_THINKING_TYPE,
    label_filter::Union{Nothing,String} = nothing,
)
    dataset = load_dataset(dataset_path)
    results_obj = load_results(results_path; dataset_path=dataset_path, provider=provider, model=model, input_fields=input_fields)
    seen = successful_result_keys(results_obj)

    processed_count = 0
    skipped_count = 0
    label_skipped_count = 0
    target_record_count = 0
    missing_input_count = 0
    normalized_label_filter = isnothing(label_filter) ? nothing :
        lowercase(strip(String(label_filter)))

    for record in dataset
        record_label = lowercase(strip(String(get(record, "label", ""))))
        if isnothing(normalized_label_filter) || record_label == normalized_label_filter
            target_record_count += 1
        end
    end

    target_record_count > 0 || throw(ArgumentError(
        "No dataset records matched label_filter=$(repr(label_filter)). " *
        "If you are evaluating newly compatible records, first finish " *
        "rewrite_task_specific_paraphrases() so they are present in the rewritten dataset.",
    ))

    println("Dataset records: ", length(dataset))
    println("Label filter: ", isnothing(label_filter) ? "none" : label_filter)
    println("Matching records: ", target_record_count)
    println("Input fields per record: ", length(input_fields))

    for record in dataset
        haskey(record, "LTL") || continue

        record_label = lowercase(strip(String(get(record, "label", ""))))
        if !isnothing(normalized_label_filter) && record_label != normalized_label_filter
            label_skipped_count += 1
            continue
        end

        for field in input_fields
            if !haskey(record, field)
                missing_input_count += 1
                @warn "Skipping missing rewritten input field" record_id=get(record, "id", nothing) field=field
                continue
            end
            source_text = strip(String(record[field]))
            isempty(source_text) && continue

            key = result_key(get(record, "id", nothing), field)
            if key in seen
                skipped_count += 1
                println("Skipping record ID ", get(record, "id", "?"), " field `", field, "` because it already has a successful result.")
                continue
            end

            println("Processing record ID ", get(record, "id", "?"), " field `", field, "`...")

            try
                predicted_ltl = request_translation(
                source_text;
                provider=provider,
                model=model,
                api_url=api_url,
                api_key_env=api_key_env,
                include_temperature=include_temperature,
                temperature=temperature,
                anthropic_version=anthropic_version,
                max_tokens=max_tokens,
                ap_mapping=get(record, "task_specific_ap_mapping", nothing),
                reasoning_effort=reasoning_effort,
                deepseek_thinking_type=deepseek_thinking_type,
)
                original_ltl = String(record["LTL"])
                equivalent = are_equivalent(predicted_ltl, original_ltl)

                entry = make_result_entry(
                    record,
                    field,
                    source_text,
                    predicted_ltl,
                    equivalent,
                    "ok";
                    translation_provider=provider,
                    translation_model=model,
                )
                append_result!(results_obj, entry)
                save_results(results_obj, results_path)
                push!(seen, key)
                processed_count += 1
                println("Saved result for record ID ", get(record, "id", "?"), " field `", field, "`. Equivalent: ", equivalent)
            catch err
                entry = make_result_entry(
                    record,
                    field,
                    source_text,
                    nothing,
                    nothing,
                    "error";
                    translation_provider=provider,
                    translation_model=model,
                    error_message=sprint(showerror, err),
                )
                append_result!(results_obj, entry)
                save_results(results_obj, results_path)
                push!(seen, key)
                processed_count += 1
                println("Saved error result for record ID ", get(record, "id", "?"), " field `", field, "`.")
            end
        end
    end

    println("Processed $(processed_count) evaluations.")
    println("Skipped $(skipped_count) evaluations already present in results.")
    println("Skipped $(label_skipped_count) records outside the label filter.")
    println("Missing rewritten input fields: $(missing_input_count).")
    println("Results path: $(results_path)")
end

function summarize_results(results_path::String = DEFAULT_RESULTS_PATH)
    results_obj = load_results(results_path; dataset_path=DEFAULT_DATASET_PATH, provider=DEFAULT_PROVIDER, model=DEFAULT_MODEL, input_fields=DEFAULT_INPUT_FIELDS)
    entries = results_obj["results"]

    total = length(entries)
    ok_count = 0
    error_count = 0
    equivalent_count = 0

    by_field = OrderedDict{String,OrderedDict{String,Int}}()

    for entry in entries
        field = haskey(entry, "input_field") ? String(entry["input_field"]) : String(entry[:input_field])
        status = haskey(entry, "status") ? String(entry["status"]) : String(entry[:status])
        equivalent = haskey(entry, "equivalent") ? entry["equivalent"] : get(entry, :equivalent, nothing)

        if !haskey(by_field, field)
            by_field[field] = OrderedDict("total" => 0, "ok" => 0, "error" => 0, "equivalent" => 0)
        end

        by_field[field]["total"] += 1
        if status == "ok"
            ok_count += 1
            by_field[field]["ok"] += 1
            if equivalent === true
                equivalent_count += 1
                by_field[field]["equivalent"] += 1
            end
        else
            error_count += 1
            by_field[field]["error"] += 1
        end
    end

    println("Results path: ", results_path)
    println("Prompt setting: ", haskey(results_obj, "prompt_setting") ? results_obj["prompt_setting"] : "unknown")
    println("Total evaluated: ", total)
    println("Successful translations: ", ok_count)
    println("Errors: ", error_count)
    println("Equivalent among successful: ", equivalent_count)
    println()

    for (field, stats) in by_field
        println("Field: ", field)
        println("  total: ", stats["total"])
        println("  ok: ", stats["ok"])
        println("  error: ", stats["error"])
        println("  equivalent: ", stats["equivalent"])
        println()
    end
end

function task_specific_success_rate(results_path::String = DEFAULT_RESULTS_PATH)
    results_obj = load_results(
        results_path;
        dataset_path=DEFAULT_DATASET_PATH,
        provider=DEFAULT_PROVIDER,
        model=DEFAULT_MODEL,
        input_fields=DEFAULT_INPUT_FIELDS,
    )

    latest = latest_results_by_key(results_obj)
    total = length(latest)
    total == 0 && begin
        println("No results found in: ", results_path)
        return 0.0
    end

    successful = 0
    equivalent = 0

    for entry in values(latest)
        status = haskey(entry, "status") ? String(entry["status"]) : ""
        if status == "ok"
            successful += 1
            if haskey(entry, "equivalent") && entry["equivalent"] === true
                equivalent += 1
            end
        end
    end

    success_rate = equivalent / total

    println("Results path: ", results_path)
    println("Latest unique evaluations: ", total)
    println("Successful translations: ", successful)
    println("Equivalent translations: ", equivalent)
    println("Success rate: ", round(success_rate * 100; digits=2), "%")

    return success_rate
end

function main()
    println("Loaded taskSpescific_evaluate.jl")
    println("Edit the configuration constants at the top of this file, then run:")
    println("  evaluate_model_on_benchmark()")
    println("or:")
    println("  summarize_results()")
    println("Temperature pilot: julia --project=. TaskTemperatureSweep.jl --run")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
