module VLTLRuntime

# Isolated modules: do not redefine the task-specific collector or its AP set.
include(joinpath(@__DIR__, "FullBenchmarkRuntime.jl"))
include(joinpath(@__DIR__, "ClusterConjunction.jl"))
const Shared = FullBenchmarkRuntime
const VLTL_DATASET = joinpath(@__DIR__, "Benchmark", "VLTL.json")
const VLTL_RESULTS = joinpath(@__DIR__, "Results", "deepseek-v4-pro-VLTL-T0.7-R10.json")
const VLTL_PROMPT_VERSION = 1

function benchmark_propositions(records)
    aps = unique([String(m.match) for r in records for m in eachmatch(r"\bprop_[1-9][0-9]*\b", String(r["NL"]))])
    isempty(aps) && error("No prop_N propositions found in the benchmark NL inputs")
    return sort!(aps; by=a->parse(Int,split(a,'_')[2]))
end

function build_vltl_prompt(sentence::AbstractString; atomic_propositions)
    aps = String.(atomic_propositions)
    isempty(aps) && error("Atomic proposition set is empty")
    names = "AP = {" * join(aps, ", ") * "}"
    # Use six generic demonstrations with the same temporal structures as the
    # task-specific prompt. No example is selected from a benchmark reference.
    a,b,c,d = [aps[min(i,length(aps))] for i in 1:4]
    return join([
        "You translate one natural-language requirement into a Linear Temporal Logic (LTL) formula.",
        "", "Task atomic propositions:", names,
        "Each prop_N is a Boolean atomic proposition at the current time.",
        "Use only the exact proposition names listed in AP. Preserve underscores and numeric suffixes.",
        "These names are propositions themselves, not aliases for colors or physical regions.",
        "", "Semantics:", "",
        "- Use standard LTL over infinite traces, evaluated from the current time.",
        "- `U` is strong Until: `(a U b)` requires `b` to eventually hold, while `a` holds at every preceding time.",
        "- “Whenever,” “every time,” and “globally” introduce a recurring global requirement.",
        "- “From then on” means that the consequent must hold globally from the triggering time.",
        "- Do not add a global operator unless the requirement expresses a global or recurring condition.",
        "", "Allowed syntax:", "",
        "- Atomic propositions: " * join(["`$ap`" for ap in aps], ", ") * ".",
        "- Negation: `!`; conjunction: `&`; disjunction: `|`; implication: `->`; equivalence: `<->`.",
        "- Next: `X`; Eventually: `F`; Globally: `G`; strong Until: `U`; Release: `R`.",
        "", "Translation rules:", "",
        "- Preserve the complete logical and temporal meaning.",
        "- Use the exact proposition names listed in AP.",
        "- Do not introduce new propositions.",
        "- Do not introduce domain assumptions, such as mutual exclusion, reachability, or ordering, unless explicitly stated.",
        "- Do not use every available proposition unless the requirement needs it.",
        "- Use parentheses to make operator scope unambiguous.",
        "- Before answering, silently check temporal scope, ordering, eventuality, and whether each Until obligation is strong.",
        "", "Output rules:", "",
        "- Output exactly one LTL formula.",
        "- Do not output explanations, labels, comments, or code fences.",
        "", "Examples:", "",
        "Natural-language requirement:", "Whenever $a holds, $b must hold from then on.",
        "", "LTL:", "G(($a) -> G($b))", "",
        "Natural-language requirement:", "Whenever $a holds, $c must hold until $b holds.",
        "", "LTL:", "G(($a) -> (($c) U ($b)))", "",
        "Natural-language requirement:", "Eventually $a holds, and from that time eventually $b holds.",
        "", "LTL:", "F(($a) & F($b))", "",
        "Natural-language requirement:", "$d must never hold.",
        "", "LTL:", "G(!($d))", "",
        "Natural-language requirement:", "$a holds until $d holds.",
        "", "LTL:", "($a) U ($d)", "",
        "Natural-language requirement:", "At the next step, $c holds.",
        "", "LTL:", "X($c)", "",
        "Target input:", "", "Available atomic proposition set:", names, "",
        "Natural-language requirement:", "<requirement>", String(sentence), "</requirement>", "", "LTL:",
    ], "\n")
end

function vltl_plan(;dataset_path=VLTL_DATASET, R=10, temperature=0.7,
    model=Shared.DEFAULT_MODEL, max_tokens=Shared.DEFAULT_MAX_TOKENS,
    api_url=Shared.DEFAULT_API_URL, api_key_env=Shared.DEFAULT_API_KEY_ENV)
    R isa Integer && !(R isa Bool) && R>0 || throw(ArgumentError("R must be a positive integer"))
    isfinite(temperature) && 0<=temperature<=2 || throw(ArgumentError("Invalid temperature"))
    max_tokens isa Integer && max_tokens>0 || throw(ArgumentError("Invalid max_tokens"))
    records = Shared.load_dataset(String(dataset_path))
    isempty(records) && error("Empty VLTL benchmark")
    all(r->haskey(r,"NL") && haskey(r,"LTL") && r["NL"] isa AbstractString &&
        r["LTL"] isa AbstractString && !isempty(strip(r["NL"])),records) || error("Expected NL,LTL records")
    aps = benchmark_propositions(records)
    prompts = Any[]
    for (i,r) in enumerate(records)
        reference = String(r["LTL"])
        isempty(setdiff(Set(conjunction_atoms(parse_ltl(reference))),Set(aps))) ||
            error("Reference $i uses propositions absent from the benchmark NL alphabet")
        push!(prompts, OrderedDict{String,Any}(
            "prompt_index"=>i, "record_id"=>i, "input_field"=>"NL",
            "source_text"=>String(r["NL"]), "ground_truth_ltl"=>reference,
            "translation_prompt"=>build_vltl_prompt(r["NL"];atomic_propositions=aps),
            "samples"=>Any[], "errors"=>Any[]))
    end
    config=OrderedDict{String,Any}(
        "benchmark"=>"VLTL", "dataset_sha256"=>bytes2hex(sha256(read(dataset_path))),
        "provider"=>"deepseek", "model"=>String(model), "temperature"=>Float64(temperature),
        "K"=>Int(R), "repetitions"=>Int(R), "input_fields"=>["NL"],
        "max_tokens"=>Int(max_tokens), "api_url"=>String(api_url), "api_key_env"=>String(api_key_env),
        "thinking_type"=>"disabled", "prompt_version"=>VLTL_PROMPT_VERSION,
        "system_prompt"=>Shared.TASK_TRANSLATION_SYSTEM_PROMPT, "atomic_propositions"=>aps,
        "prompt_policy"=>"generic Boolean prop_N; exact allowed set; no reference in request",
        "adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "response_policy"=>"retain_all_returned_outputs_without_logical_filtering")
    return OrderedDict{String,Any}("schema_version"=>1,"stage"=>"LLM_collection_only",
        "dataset_path"=>abspath(dataset_path), "config"=>config, "record_count"=>length(prompts),
        "prompt_count"=>length(prompts), "planned_queries"=>length(prompts)*R,
        "created_at_utc"=>Shared.utc_timestamp(), "updated_at_utc"=>Shared.utc_timestamp(),
        "run_state"=>"prepared", "prompts"=>prompts)
end

function prepare_collection(path;kwargs...)
    expected=vltl_plan(;kwargs...)
    if !isfile(path)
        Shared.checkpoint!(expected,path)
        return expected
    end
    saved=Shared.mutable_json(JSON3.read(read(path,String)))
    saved["config"]==expected["config"] || error("VLTL configuration/dataset changed; use a new output path")
    saved["stage"]==expected["stage"] && saved["prompt_count"]==expected["prompt_count"] &&
        saved["planned_queries"]==expected["planned_queries"] || error("Invalid VLTL checkpoint")
    length(saved["prompts"])==length(expected["prompts"]) || error("Prompt count changed")
    for (p,q) in zip(saved["prompts"],expected["prompts"])
        for field in keys(q)
            field in ("samples","errors") && continue
            p[field]==q[field] || error("Saved VLTL prompt differs at $field")
        end
        length(p["samples"])<=expected["config"]["K"] || error("Too many saved samples")
        for (j,s) in enumerate(p["samples"])
            s["query_index"]==s["repetition"]==j && s["status"]=="returned" || error("Invalid query sequence")
            s["formula"] isa AbstractString && s["raw_output"] isa AbstractString || error("Invalid saved response")
        end
    end
    return saved
end

function prepare_vltl(;results_path=VLTL_RESULTS,kwargs...)
    path=abspath(results_path)
    return Shared.with_collection_lock(path) do
        prepare_collection(path;kwargs...)
    end
end

vltl_status(;results_path=VLTL_RESULTS)=Shared.collection_summary(Shared.mutable_json(JSON3.read(read(results_path,String))))

function request_body(prompt,config)
    # This pure function is also checked offline: no color-specific builder.
    return Dict("model"=>config["model"], "temperature"=>config["temperature"],
        "max_tokens"=>config["max_tokens"], "stream"=>false,"thinking"=>Dict("type"=>"disabled"),
        "messages"=>[Dict("role"=>"system","content"=>config["system_prompt"]),
                     Dict("role"=>"user","content"=>String(prompt))])
end

function request_vltl(prompt,config)
    api_key=Shared.get_api_key(config["api_key_env"])
    body=request_body(prompt,config)
    response=Shared.post_json(config["api_url"],
        ["Authorization"=>"Bearer $api_key","Content-Type"=>"application/json"],JSON3.write(body))
    response.status==200 || error("DeepSeek API request failed with status $(response.status): $(String(response.body))")
    parsed=JSON3.read(String(response.body))
    return Shared.translation_response(parsed,Shared.extract_deepseek_output_text(parsed);return_metadata=true)
end

function run_vltl(;results_path=VLTL_RESULTS,workers=2,request_fn=request_vltl,verbose=true,kwargs...)
    workers isa Integer && workers>0 || throw(ArgumentError("workers must be positive"))
    path=abspath(results_path)
    return Shared.with_collection_lock(path) do
        data=prepare_collection(path;kwargs...)
        Shared.collection_summary(data)["complete"] && return data["summary"]
        config=data["config"]
        request_fn===request_vltl && Shared.get_api_key(config["api_key_env"])
        jobs=Channel{Int}(data["prompt_count"])
        for (i,p) in enumerate(data["prompts"])
            length(p["samples"])<config["K"] && put!(jobs,i)
        end
        close(jobs)
        mutex=ReentrantLock();stopped=Ref(false);failure=Ref{Any}(nothing);tasks=Task[]
        data["run_state"]="running";Shared.checkpoint!(data,path)
        try
            @sync for _ in 1:min(workers,data["prompt_count"])
                task=@async try
                    for i in jobs
                        stopped[] && break
                        p=data["prompts"][i]
                        for j in (length(p["samples"])+1):config["K"]
                            stopped[] && break
                            start=time()
                            response=try
                                request_fn(p["translation_prompt"],config)
                            catch err
                                err isa InterruptException && rethrow()
                                lock(mutex) do
                                    push!(p["errors"],Dict("query_index"=>j,"error"=>sprint(showerror,err),"at_utc"=>Shared.utc_timestamp()))
                                    failure[]===nothing && (failure[]=sprint(showerror,err))
                                    stopped[]=true
                                end
                                break
                            end
                            lock(mutex) do
                                push!(p["samples"],OrderedDict{String,Any}(
                                    "query_index"=>j,"repetition"=>j,"status"=>"returned",
                                    "formula"=>String(response.formula),"raw_output"=>String(response.raw_output),
                                    "usage"=>Shared.mutable_json(response.usage),"returned_model"=>response.returned_model,
                                    "finish_reason"=>response.finish_reason,"elapsed_seconds"=>time()-start,
                                    "received_at_utc"=>Shared.utc_timestamp()))
                                s=Shared.checkpoint!(data,path) # save every returned response
                                verbose && println("VLTL saved $(s["saved_responses"])/$(s["planned_queries"]) responses")
                                flush(stdout)
                            end
                        end
                    end
                catch
                    stopped[]=true;rethrow()
                end
                push!(tasks,task)
            end
        catch
            stopped[]=true
            for t in tasks
                try wait(t) catch end
            end
            rethrow()
        finally
            data["run_state"]=failure[]!==nothing ? "failed" : Shared.collection_summary(data)["complete"] ? "complete" : "interrupted"
            Shared.checkpoint!(data,path)
        end
        failure[]===nothing || error("VLTL collection stopped; saved responses preserved. Rerun run_vltl() to resume. First error: $(failure[])")
        return data["summary"]
    end
end

function analyze_vltl(;source_path=VLTL_RESULTS,results_directory=splitext(source_path)[1]*"-conjunctions",max_k=10)
    data=conjunction_read(source_path);config=data["config"]
    config["benchmark"]=="VLTL" || error("Expected a VLTL collection")
    max_k isa Integer && 1<=max_k<=config["K"] || error("Require 1 <= max_k <= repetitions")
    data["summary"]["complete"] || error("Finish VLTL collection before semantic analysis")
    aps=Set(String.(config["atomic_propositions"]))
    length(data["prompts"])==data["prompt_count"] || error("Prompt count mismatch")
    for (i,p) in enumerate(data["prompts"])
        p["prompt_index"]==i || error("Invalid prompt ordering")
        length(p["samples"])==config["K"] || error("Incomplete query set at prompt $i")
        [s["query_index"] for s in p["samples"]]==collect(1:config["K"]) || error("Invalid query sequence")
        all(s->s["status"]=="returned",p["samples"]) || error("Non-returned sample")
        isempty(setdiff(Set(conjunction_atoms(parse_ltl(p["ground_truth_ltl"]))),aps)) || error("Invalid reference alphabet")
    end
    mkpath(results_directory)
    identity=Dict("source_sha256"=>bytes2hex(sha256(read(source_path))),
        "adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "engine_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"ClusterConjunction.jl")))) ,
        "max_k"=>max_k,"candidate_filter"=>"parse_then_allowed_APs_then_SAT_before_clustering",
        "ranking"=>"frequency_descending_then_earliest_query","missing_policy"=>"saturate",
        "denominator_policy"=>"all_prompts","allowed_propositions"=>sort!(collect(aps)))
    return Shared.with_collection_lock(joinpath(results_directory,"summary.json")) do
        output=joinpath(results_directory,"summary.json")
        completed=Any[]
        if isfile(output)
            previous=conjunction_read(output)
            previous["analysis_identity"]==identity || error("Analysis/source changed; use a new results directory")
            append!(completed,previous["prompt_results"])
        end
        report=OrderedDict{String,Any}("analysis_identity"=>identity,"source_path"=>abspath(source_path),
            "benchmark"=>"VLTL","model"=>config["model"],"temperature"=>config["temperature"],
            "repetitions"=>config["K"],"prompt_count"=>data["prompt_count"],
            "samples"=>data["planned_queries"],"prompt_results"=>completed,"complete"=>false)
        for p in data["prompts"][(length(completed)+1):end]
            reference=String(p["ground_truth_ltl"])
            samples=[Dict("candidate_formula"=>s["formula"],"repetition"=>s["query_index"],
                "finish_reason"=>s["finish_reason"]) for s in p["samples"]]
            excluded=Any[]
            clusters=rank_conjunction_clusters(samples,reference;excluded_candidates=excluded,allowed_propositions=aps)
            push!(completed,OrderedDict{String,Any}("prompt_index"=>p["prompt_index"],"record_id"=>p["record_id"],
                "input_field"=>"NL","source_text"=>p["source_text"],"reference_ltl"=>reference,
                "samples"=>length(samples),"satisfiable_samples"=>sum(c["frequency"] for c in clusters;init=0),
                "excluded_candidates"=>excluded,"ranked_clusters"=>clusters,
                "conjunctions"=>evaluate_cluster_conjunctions(clusters,reference;max_k=max_k)))
            report["completed_prompts"]=length(completed)
            conjunction_save(output,report)
            println("VLTL analyzed $(length(completed))/$(data["prompt_count"]) prompts");flush(stdout)
        end
        report["complete"]=length(completed)==data["prompt_count"]
        report["summary"]=aggregate_all_prompt_conjunctions(completed;max_k=max_k)
        excluded=[c for p in completed for c in p["excluded_candidates"]]
        report["exclusion_counts"]=Dict(name=>count(c->c["reason"]==name,excluded)
            for name in ("invalid_ltl","disallowed_propositions","unsatisfiable"))
        report["completed_prompts"]=length(completed)
        conjunction_save(output,report)
        open(joinpath(results_directory,"comparison.md"),"w") do io
            println(io,"# VLTL: SAT-filtered semantic inclusion\n\n42 source pairs; scores use the supplied references without modification.\n")
            println(io,"| k | Logical good | Not guaranteed | Abstain | Unsatisfiable |\n|---|---:|---:|---:|---:|")
            for r in report["summary"]
                cells=[@sprintf("%.2f%%",100*r["rates"][name]) for name in ("logical_good","not_guaranteed","abstain","unsatisfiable")]
                println(io,"| $(r["k"]) | $(join(cells," | ")) |")
            end
        end
        return report
    end
end

end # module
