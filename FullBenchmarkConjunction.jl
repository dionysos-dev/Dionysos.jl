if !isdefined(@__MODULE__, :FullBenchmarkConjunctionRuntime)
@eval module FullBenchmarkConjunctionRuntime

include(joinpath(@__DIR__, "ClusterConjunction.jl"))

const FULL_CONJUNCTION_SOURCE = joinpath(@__DIR__, "Results", "deepseek-v4-pro-task-specific-full-T0.7-K7.json")

function source_prompts(source)
    data = conjunction_read(source)
    get(data,"stage",nothing)=="LLM_collection_only" || error("Expected a full-benchmark collection JSON")
    config = data["config"]
    config["prompt_version"]==2 || error("Expected task-specific prompt version 2")
    Set(config["atomic_propositions"])==CONJUNCTION_ALLOWED_PROPOSITIONS || error("Atomic proposition set differs")
    prompts = data["prompts"]
    isempty(prompts) && error("Source has no prompts")
    length(prompts)==data["prompt_count"] || error("Prompt count mismatch")
    K = Int(config["K"])
    data["planned_queries"]==K*length(prompts) || error("Planned query count mismatch")
    for (i,p) in enumerate(prompts)
        p["prompt_index"]==i || error("Prompt indices must be contiguous and ordered")
        samples=p["samples"]
        length(samples)==K || error("Incomplete prompt $i: expected $K samples")
        [s["query_index"] for s in samples]==collect(1:K) || error("Invalid query sequence at prompt $i")
        all(s -> s["status"]=="returned",samples) || error("Non-returned sample at prompt $i")
        reference=p["ground_truth_ltl"]
        reference isa AbstractString || error("Missing reference at prompt $i")
        ast=parse_ltl(reference)
        isempty(setdiff(Set(conjunction_atoms(ast)),CONJUNCTION_ALLOWED_PROPOSITIONS)) ||
            error("Reference has disallowed propositions at prompt $i")
    end
    return data
end

function analyze_prompt(p, max_k)
    # Adapt acquisition fields only; the frozen development filtering and logic
    # functions below are reused unchanged. No source outputs are overwritten.
    samples = [OrderedDict{String,Any}("candidate_formula"=>s["formula"],
        "repetition"=>s["query_index"],"finish_reason"=>s["finish_reason"])
        for s in p["samples"]]
    excluded=Any[]
    reference=String(p["ground_truth_ltl"])
    clusters=rank_conjunction_clusters(samples,reference;excluded_candidates=excluded)
    return OrderedDict{String,Any}(
        "prompt_index"=>p["prompt_index"],"record_id"=>p["record_id"],
        "input_field"=>p["input_field"],"source_text"=>p["source_text"],
        "development_prompt_index"=>get(p,"development_prompt_index",nothing),
        "reference_ltl"=>reference,"samples"=>length(samples),
        "satisfiable_samples"=>sum(c["frequency"] for c in clusters;init=0),
        "excluded_candidates"=>excluded,"ranked_clusters"=>clusters,
        "conjunctions"=>evaluate_cluster_conjunctions(clusters,reference;max_k=max_k))
end

function write_report(directory, report)
    metrics=["logical_good","not_guaranteed","abstain","unsatisfiable","top1_rescue","semantic_rescue"]
    open(joinpath(directory,"comparison.md"),"w") do io
        println(io,"# Full-benchmark SAT-filtered conjunction curves\n")
        println(io,"$(report["model"]), T=$(report["temperature"]), $(report["prompt_count"]) prompts, $(report["samples"]) saved responses. Existing outputs only; no model/API calls.\n")
        println(io,"Parse → allowed APs → individual SAT → discard failures → semantic equivalence clusters → rank by frequency (earliest query breaks ties) → conjoin min(k,m) representatives. All rates use all prompts; zero clusters means abstention.\n")
        println(io,"| k | P_good | P_not-guaranteed | P_abstain | P_unsat | Top-1 rescue | Semantic rescue | Representatives used |")
        println(io,"|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in report["summary"]
            values=[@sprintf("%.2f%%",100*row["rates"][m]) for m in metrics]
            println(io,"| $(row["k"]) | $(join(values," | ")) | $(row["included_representatives"]) |")
        end
        println(io,"\nGood means SAT(Φ_k) and Φ_k implies the reference, over unconstrained infinite traces. Unsatisfiable conjunctions do not count as good.\n")
        println(io,"Top-1 rescue means Good(Φ_k) and the top representative is not sufficient. Semantic rescue means Good(Φ_k) and no included representative is equivalent to the reference. Sufficient means SAT(candidate) and candidate implies reference. Representatives used is the sum of min(k,m) across all prompts.\n")
        println(io,"Excluded candidates: $(report["exclusion_counts"]). Retained samples: $(first(report["summary"])["available_candidates"]); semantic clusters summed over prompts: $(first(report["summary"])["available_semantic_clusters"]).\n")
        println(io,"These are logical metrics; controller synthesis is not part of this analysis. All counts, cluster members, formulas, per-prompt outcomes, source SHA-256, and policy details are saved in summary.json.")
    end
end

function analyze_full_benchmark_conjunctions(; source_path=FULL_CONJUNCTION_SOURCE,
    results_directory=splitext(source_path)[1]*"-conjunctions", max_k=7, workers=4,
    checkpoint_every=50)
    max_k isa Integer && max_k>0 || throw(ArgumentError("max_k must be positive"))
    workers isa Integer && workers>0 || throw(ArgumentError("workers must be positive"))
    checkpoint_every isa Integer && checkpoint_every>0 || throw(ArgumentError("checkpoint_every must be positive"))
    Sys.which("ltlfilt")===nothing && error("Spot's ltlfilt is required")
    data=source_prompts(source_path)
    prompts=data["prompts"]
    identity=OrderedDict{String,Any}(
        "schema_version"=>1,"source_sha256"=>bytes2hex(sha256(read(source_path))),
        "adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "engine_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"ClusterConjunction.jl")))),
        "max_k"=>max_k,"candidate_filter"=>"parse_then_allowed_APs_then_SAT_before_clustering",
        "ranking"=>"frequency_descending_then_earliest_query",
        "missing_policy"=>"saturate","denominator_policy"=>"all_prompts",
        "rescue_definition"=>"no_individual_equivalent")
    mkpath(results_directory)
    lock_path=joinpath(results_directory,".analysis.lock")
    try
        mkdir(lock_path)
    catch
        isdir(lock_path) || rethrow()
        error("Analysis lock exists: $lock_path. Check for an active analysis before removing a stale lock.")
    end
    try
        output=joinpath(results_directory,"summary.json")
        completed=Any[nothing for _ in prompts]
        if isfile(output)
            saved=conjunction_read(output)
            saved["analysis_identity"]==identity || error("Analysis/source changed; use a new results_directory")
            for p in saved["prompt_results"]
                i=Int(p["prompt_index"])
                1<=i<=length(prompts) && completed[i]===nothing || error("Invalid checkpoint prompt index")
                completed[i]=p
            end
        end
        report=OrderedDict{String,Any}("analysis_identity"=>identity,
            "source_path"=>abspath(source_path),"model"=>data["config"]["model"],
            "temperature"=>data["config"]["temperature"],"repetitions"=>data["config"]["K"],
            "prompt_count"=>length(prompts),"samples"=>sum(length(p["samples"]) for p in prompts))
        function checkpoint()
            report["updated_at_utc"]=string(now(UTC))*"Z"
            report["prompt_results"]=Any[p for p in completed if p!==nothing]
            report["completed_prompts"]=length(report["prompt_results"])
            report["complete"]=report["completed_prompts"]==length(prompts)
            if report["complete"]
                report["summary"]=aggregate_all_prompt_conjunctions(report["prompt_results"];max_k=max_k)
                excluded=[c for p in report["prompt_results"] for c in p["excluded_candidates"]]
                report["exclusion_counts"]=OrderedDict{String,Any}(name=>count(c->c["reason"]==name,excluded)
                    for name in ("invalid_ltl","disallowed_propositions","unsatisfiable"))
            end
            conjunction_save(output,report)
            println("Saved semantic analysis: $(report["completed_prompts"])/$(length(prompts)) prompts")
            flush(stdout)
        end
        pending=[i for i in eachindex(prompts) if completed[i]===nothing]
        jobs=Channel{Int}(length(pending))
        for i in pending
            put!(jobs,i)
        end
        close(jobs)
        state_lock=ReentrantLock()
        stopped=Ref(false)
        failure=Ref{Any}(nothing)
        dirty=Ref(0)
        tasks=Task[]
        try
            @sync for _ in 1:min(workers,length(pending))
                task=@async try
                    for i in jobs
                        stopped[] && break
                        result=try
                            analyze_prompt(prompts[i],max_k)
                        catch err
                            stopped[]=true
                            failure[]=(prompt_index=i,error=sprint(showerror,err))
                            rethrow()
                        end
                        lock(state_lock) do
                            completed[i]=result
                            dirty[]+=1
                            if dirty[]>=checkpoint_every
                                checkpoint()
                                dirty[]=0
                            end
                        end
                    end
                catch
                    stopped[]=true
                    rethrow()
                end
                push!(tasks,task)
            end
        catch
            stopped[]=true
            for task in tasks
                try
                    wait(task)
                catch
                end
            end
            rethrow()
        finally
            failure[]===nothing || (report["last_error"]=failure[])
            checkpoint()
        end
        write_report(results_directory,report)
        println("Table: ",joinpath(results_directory,"comparison.md"))
        return report
    finally
        rm(lock_path;recursive=true)
    end
end

end # module
end # guard

analyze_full_benchmark_conjunctions(;kwargs...)=FullBenchmarkConjunctionRuntime.analyze_full_benchmark_conjunctions(;kwargs...)

if abspath(PROGRAM_FILE)==@__FILE__
    isempty(ARGS) || ARGS==["--run"] || error("Use --run (saved outputs only; no model calls)")
    analyze_full_benchmark_conjunctions()
end
