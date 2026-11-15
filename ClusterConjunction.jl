using JSON3
using OrderedCollections: OrderedDict
using Dates
using SHA
using Printf

if !isdefined(@__MODULE__, :parse_ltl)
    include(joinpath(@__DIR__, "LTL.jl"))
end

const CONJUNCTION_TEMPERATURE = 0.7
const CONJUNCTION_SOURCE = joinpath(@__DIR__, "Results", "task-specific-temperature-prompt-v2")
const CONJUNCTION_ALLOWED_PROPOSITIONS = Set(["blue","brown","green","purple","yellow"])

conjunction_mutable(x::AbstractDict) = OrderedDict{String,Any}(String(k)=>conjunction_mutable(v) for (k,v) in pairs(x))
conjunction_mutable(x::AbstractVector) = Any[conjunction_mutable(v) for v in x]
conjunction_mutable(x) = x
conjunction_read(path) = conjunction_mutable(JSON3.read(read(path, String)))

function conjunction_save(path, data)
    mkpath(dirname(abspath(path)))
    temporary, io = mktemp(dirname(abspath(path)))
    try
        JSON3.pretty(io, data)
        write(io, '\n')
        close(io)
        mv(temporary, path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

"""Use Spot on infinite traces. Parse/tool errors and timeouts never become UNSAT."""
function conjunction_spot_match(formula, filter; timeout_seconds=120.0)
    executable = Sys.which("ltlfilt")
    isnothing(executable) && throw(ArgumentError("Spot's ltlfilt is required"))
    out, err = IOBuffer(), IOBuffer()
    process = run(pipeline(ignorestatus(Cmd([executable, "-f", String(formula),
        String(filter), "--count"])); stdout=out, stderr=err); wait=false)
    deadline = time() + timeout_seconds
    while process_running(process)
        if time() > deadline
            kill(process)
            wait(process)
            error("Spot timed out after $timeout_seconds seconds for $filter")
        end
        sleep(0.01)
    end
    wait(process)
    process.exitcode in (0,1) || error("Spot failed: $(strip(String(take!(err))))")
    return process.exitcode == 0
end

const CONJUNCTION_SAT_CACHE = Dict{String,Bool}()
const CONJUNCTION_IMPLIES_CACHE = Dict{Tuple{String,String},Bool}()
const CONJUNCTION_EQUIVALENCE_CACHE = Dict{Tuple{String,String},Bool}()
conjunction_sat(f) = get!(() -> !conjunction_spot_match(f, "--equivalent-to=0"), CONJUNCTION_SAT_CACHE, String(f))
conjunction_implies(f, reference) = get!(() -> conjunction_spot_match(f, "--imply=$reference"),
    CONJUNCTION_IMPLIES_CACHE, (String(f),String(reference)))
function conjunction_equivalent(a,b)
    a == b && return true
    key = a < b ? (String(a),String(b)) : (String(b),String(a))
    return get!(() -> conjunction_spot_match(a,"--equivalent-to=$b"), CONJUNCTION_EQUIVALENCE_CACHE, key)
end

conjunction_atoms(node::LTLAtom) = node.is_constant ? String[] :
    [startswith(node.name,'"') ? String(JSON3.read(node.name)) : node.name]
conjunction_atoms(node::LTLUnary) = conjunction_atoms(node.child)
conjunction_atoms(node::LTLBinary) = vcat(conjunction_atoms(node.left),conjunction_atoms(node.right))

"""Parse, validate APs, then check SAT for every candidate before clustering."""
function filter_conjunction_candidates(samples; excluded_candidates=Any[],
    allowed_propositions=CONJUNCTION_ALLOWED_PROPOSITIONS)
    eligible = Any[]
    for sample in samples
        candidate = get(sample,"candidate_formula",nothing)
        ast = try
            candidate isa AbstractString || throw(ArgumentError("Missing formula string"))
            parse_ltl(candidate)
        catch err
            err isa ArgumentError || rethrow()
            push!(excluded_candidates,OrderedDict{String,Any}("repetition"=>sample["repetition"],
                "reason"=>"invalid_ltl","candidate_formula"=>candidate,
                "validation_error"=>sprint(showerror,err)))
            continue
        end
        formula = String(candidate)
        unknown = setdiff(Set(conjunction_atoms(ast)),Set(allowed_propositions))
        if !isempty(unknown)
            push!(excluded_candidates,OrderedDict{String,Any}("repetition"=>sample["repetition"],
                "reason"=>"disallowed_propositions","candidate_formula"=>formula,
                "disallowed_propositions"=>sort!(collect(unknown))))
            continue
        end
        # Retain the acquisition layer's rejection of incomplete/API responses.
        if !get(sample,"valid",true) || get(sample,"finish_reason",nothing)=="length"
            push!(excluded_candidates,OrderedDict{String,Any}("repetition"=>sample["repetition"],
                "reason"=>"invalid_ltl","candidate_formula"=>formula))
            continue
        end
        if !conjunction_sat(formula)
            push!(excluded_candidates,OrderedDict{String,Any}("repetition"=>sample["repetition"],
                "reason"=>"unsatisfiable","candidate_formula"=>formula))
            continue
        end
        push!(eligible,sample)
    end
    return eligible
end

"""Cluster only the completely filtered candidate pool, then rank by frequency."""
function rank_conjunction_clusters(samples, reference; excluded_candidates=Any[],
    allowed_propositions=CONJUNCTION_ALLOWED_PROPOSITIONS)
    eligible = filter_conjunction_candidates(samples;excluded_candidates=excluded_candidates,
        allowed_propositions=allowed_propositions)
    clusters = Any[]
    for sample in eligible
        formula = String(sample["candidate_formula"])
        index = findfirst(c -> conjunction_equivalent(formula,c["representative"]), clusters)
        if isnothing(index)
            push!(clusters, OrderedDict{String,Any}(
                "cluster_id"=>length(clusters)+1, "representative"=>formula,
                "frequency"=>0, "sample_indices"=>Int[],
            ))
            index = length(clusters)
        end
        cluster = clusters[index]
        cluster["frequency"] += 1
        push!(cluster["sample_indices"],Int(sample["repetition"]))
    end
    sort!(clusters; by=c -> (-c["frequency"],minimum(c["sample_indices"])))
    for (rank,cluster) in enumerate(clusters)
        formula = cluster["representative"]
        sat = conjunction_sat(formula)
        implies = !sat || conjunction_implies(formula,reference)
        cluster["rank"] = rank
        cluster["satisfiable"] = sat
        cluster["implies_reference"] = implies
        cluster["individually_sufficient"] = sat && implies
        cluster["equivalent_to_reference"] = conjunction_equivalent(formula,reference)
    end
    return clusters
end

"""Evaluate exact prefixes Φ_k; no repeated or invented clusters are appended."""
function evaluate_cluster_conjunctions(clusters, reference; max_k=10)
    rows = Any[]
    max_k >= 1 || throw(ArgumentError("max_k must be positive"))
    for k in 1:min(max_k,length(clusters))
        included = clusters[1:k]
        formula = join(["($(c["representative"]))" for c in included]," & ")
        # A contradiction persists under additional conjunctions. Implication
        # also persists, although satisfiability may disappear at the next rank.
        previous = isempty(rows) ? nothing : last(rows)
        sat = ((!isnothing(previous) && !previous["satisfiable"]) ||
            any(c -> !c["satisfiable"],included)) ? false : conjunction_sat(formula)
        implies = !sat || (!isnothing(previous) && previous["implies_reference"]) ||
            any(c -> c["implies_reference"],included) || conjunction_implies(formula,reference)
        good = sat && implies
        push!(rows, OrderedDict{String,Any}(
            "k"=>k, "phi_k"=>formula, "cluster_ranks"=>collect(1:k),
            "satisfiable"=>sat, "implies_reference"=>implies,
            "logical_good"=>good, "not_guaranteed"=>sat && !implies,
            "unsatisfiable"=>!sat,
            "semantic_rescue_no_individual_sufficient"=>good && !any(c -> c["individually_sufficient"],included),
            "semantic_rescue_top_not_sufficient"=>good && !first(clusters)["individually_sufficient"],
            "semantic_rescue_no_individual_equivalent"=>good && !any(c -> c["equivalent_to_reference"],included),
        ))
    end
    return rows
end

function aggregate_cluster_conjunctions(prompt_results; max_k=10,
    missing_policy="saturate", rescue_definition="no_individual_equivalent")
    missing_policy in ("strict","saturate") || throw(ArgumentError("missing_policy must be strict or saturate"))
    rescue_key = "semantic_rescue_" * rescue_definition
    rescue_definition in ("no_individual_sufficient","top_not_sufficient","no_individual_equivalent") ||
        throw(ArgumentError("Unknown semantic-rescue definition"))
    rows = Any[]
    for k in 1:max_k
        selected = Any[]
        for prompt in prompt_results
            prefixes = prompt["conjunctions"]
            isempty(prefixes) && continue # Never substitute true for a failed translation.
            if length(prefixes) >= k
                push!(selected,prefixes[k])
            elseif missing_policy == "saturate"
                push!(selected,last(prefixes))
            end
        end
        n = length(selected)
        counts = OrderedDict{String,Any}(name=>count(r -> r[key],selected) for (name,key) in [
            ("satisfiable","satisfiable"),("implies_reference","implies_reference"),
            ("logical_good","logical_good"),("not_guaranteed","not_guaranteed"),
            ("unsatisfiable","unsatisfiable"),("semantic_rescue",rescue_key),
            ("top1_rescue","semantic_rescue_top_not_sufficient"),
            ("synergistic_rescue","semantic_rescue_no_individual_sufficient")])
        @assert counts["logical_good"]+counts["not_guaranteed"]+counts["unsatisfiable"] == n
        @assert counts["satisfiable"]+counts["unsatisfiable"] == n
        @assert counts["implies_reference"] == counts["logical_good"]+counts["unsatisfiable"]
        @assert counts["semantic_rescue"] <= counts["logical_good"]
        @assert counts["synergistic_rescue"] <= counts["top1_rescue"] <= counts["logical_good"]
        rates = OrderedDict{String,Any}(name=>(n==0 ? nothing : count/n) for (name,count) in counts)
        push!(rows, OrderedDict{String,Any}("k"=>k,"eligible_prompts"=>n,
            "total_prompts"=>length(prompt_results),"counts"=>counts,"rates"=>rates))
    end
    return rows
end

"""Unconditional saturated curves: zero retained clusters mean abstention.

Available candidates count retained samples, including duplicates. Available
semantic clusters count their equivalence classes; included representatives
count the actual min(k,m) terms used across all prompts.
"""
function aggregate_all_prompt_conjunctions(prompt_results; max_k=5,
    rescue_definition="no_individual_equivalent")
    rows = aggregate_cluster_conjunctions(prompt_results; max_k=max_k,
        missing_policy="saturate", rescue_definition=rescue_definition)
    total = length(prompt_results)
    available = sum(p["satisfiable_samples"] for p in prompt_results; init=0)
    clusters = sum(length(p["ranked_clusters"]) for p in prompt_results; init=0)
    for row in rows
        row["counts"]["abstain"] = total - row["eligible_prompts"]
        row["rates"] = OrderedDict{String,Any}(name=>(total==0 ? nothing : count/total)
            for (name,count) in row["counts"])
        row["denominator_policy"] = "all_prompts"
        row["available_candidates"] = available
        row["available_semantic_clusters"] = clusters
        row["included_representatives"] = sum(min(row["k"],length(p["ranked_clusters"])) for p in prompt_results; init=0)
        row["prompts_with_at_least_k_clusters"] = count(p -> length(p["ranked_clusters"])>=row["k"],prompt_results)
        @assert sum(row["counts"][name] for name in ("logical_good","not_guaranteed","abstain","unsatisfiable"))==total
    end
    return rows
end

function write_conjunction_tables(directory, analysis; missing_policy="saturate",
    rescue_definition="no_individual_equivalent")
    rows = aggregate_cluster_conjunctions(analysis["prompt_results"];
        max_k=analysis["max_k"],missing_policy=missing_policy,rescue_definition=rescue_definition)
    report = merge(analysis, OrderedDict{String,Any}("missing_policy"=>missing_policy,
        "rescue_definition"=>rescue_definition,"summary"=>rows))
    conjunction_save(joinpath(directory,"summary.json"),report)
    names = ["satisfiable","implies_reference","logical_good","not_guaranteed","unsatisfiable","semantic_rescue","top1_rescue","synergistic_rescue"]
    open(joinpath(directory,"summary.csv"),"w") do io
        println(io,join(vcat(["k","eligible_prompts","total_prompts"],names,[n*"_count" for n in names]),','))
        for row in rows
            values = vcat(Any[row["k"],row["eligible_prompts"],row["total_prompts"]],
                Any[isnothing(row["rates"][n]) ? "" : row["rates"][n] for n in names],
                Any[row["counts"][n] for n in names])
            println(io,join(values,','))
        end
    end
    open(joinpath(directory,"comparison.md"),"w") do io
        setting = isnothing(analysis["temperature"]) ? "provider-default temperature" : "T=$(analysis["temperature"])"
        model = get(analysis,"model",nothing)
        println(io,"Frequency-ranked cluster conjunctions — $(isnothing(model) ? "" : "$model, ")$setting\n")
        println(io,"$(length(analysis["prompt_results"])) prompts; $(analysis["samples"]) saved samples. No new LLM calls. Candidates must be valid and individually satisfiable before semantic clustering; ties use the earliest sampled member.\n")
        excluded = [c for p in analysis["prompt_results"] for c in p["excluded_candidates"]]
        println(io,"Filtered $(count(c -> c["reason"]=="invalid_ltl",excluded)) invalid candidates, $(count(c -> c["reason"]=="disallowed_propositions",excluded)) candidates with disallowed propositions, and $(count(c -> c["reason"]=="unsatisfiable",excluded)) unsatisfiable candidates before clustering.\n")
        println(io,"| k | Eligible prompts | Satisfiable | Implies reference | Logical good | Not guaranteed | Unsatisfiable | Semantic rescue | Top-1 rescue | Synergistic rescue |")
        println(io,"|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in rows
            values = [isnothing(row["rates"][n]) ? "—" : @sprintf("%.1f%%",100*row["rates"][n]) for n in names]
            println(io,"| $(row["k"]) | $(row["eligible_prompts"]) | $(join(values," | ")) |")
        end
        println(io,"\nMissing-cluster policy: $missing_policy. Strict rows estimate probabilities conditional on at least k distinct clusters; their populations change with k. Saturated rows conjoin min(k,m) clusters. Prompts with zero eligible clusters are unavailable under both policies; eligible means valid and individually satisfiable.\n")
        println(io,"Semantic rescue definition: $rescue_definition. The selected default is SAT(Φ_k) and Φ_k implies the reference, while no included cluster is individually equivalent to the reference. This includes stronger sufficient formulas, so rescue can be nonzero at k=1.\n")
        println(io,"Top-1 rescue: Good(Φ_k) and the top-ranked representative is not individually sufficient. Synergistic rescue: Good(Φ_k) and none of the included representatives is individually sufficient. Here sufficient means SAT(ψ) and ψ implies the reference. Both rescue metrics are zero at k=1. Synergistic rescue can disappear when an individually sufficient cluster is added, even if the conjunction stays logically good.\n")
        println(io,"Raw implication includes vacuous truth for unsatisfiable conjunctions. Logical good explicitly requires satisfiability. Logical good + not guaranteed + unsatisfiable = 100% for every nonempty row.\n")
        println(io,"All checks use standard LTL over unconstrained infinite traces. Every ranked representative, frequency, and Φ_k is saved in summary.json. A second table with the alternative missing-cluster policy is in summary-alternative-policy.json.")
    end
    open(joinpath(directory,"rescue-metrics.csv"),"w") do io
        println(io,"k,eligible_prompts,top1_rescue,top1_rescue_count,synergistic_rescue,synergistic_rescue_count")
        for row in rows
            values = Any[row["k"],row["eligible_prompts"],
                isnothing(row["rates"]["top1_rescue"]) ? "" : row["rates"]["top1_rescue"],row["counts"]["top1_rescue"],
                isnothing(row["rates"]["synergistic_rescue"]) ? "" : row["rates"]["synergistic_rescue"],row["counts"]["synergistic_rescue"]]
            println(io,join(values,','))
        end
    end
    alternate_policy = missing_policy=="strict" ? "saturate" : "strict"
    conjunction_save(joinpath(directory,"summary-alternative-policy.json"),OrderedDict{String,Any}(
        "missing_policy"=>alternate_policy,"rescue_definition"=>rescue_definition,
        "summary"=>aggregate_cluster_conjunctions(analysis["prompt_results"];
            max_k=analysis["max_k"],missing_policy=alternate_policy,rescue_definition=rescue_definition)))
    return report
end

"""Analyze saved version-2 samples at the selected temperature (nothing for default).

No model/API access. Prompts with no retained clusters have no conjunction and
are reported as unavailable rather than counted as an empty true formula.
"""
function analyze_cluster_conjunctions(; source_directory=CONJUNCTION_SOURCE,
    temperature=CONJUNCTION_TEMPERATURE,
    results_directory=joinpath(source_directory,"conjunction-temperature-$(isnothing(temperature) ? "default" : @sprintf("%.1f",temperature))-validated"),
    max_k=10, missing_policy="saturate",rescue_definition="no_individual_equivalent")
    manifest = conjunction_read(joinpath(source_directory,"manifest.json"))
    get(manifest,"prompt_version",0)==2 || error("Use the corrected version-2 prompt run")
    temperature in manifest["temperatures"] || error("Source run has no samples for temperature=$temperature")
    isnothing(temperature) && get(manifest,"include_temperature",true) &&
        error("Default-temperature checkpoints require include_temperature=false")
    Set(manifest["atomic_propositions"])==CONJUNCTION_ALLOWED_PROPOSITIONS || error("Source atomic proposition set differs from the task vocabulary")
    identity = OrderedDict{String,Any}("schema_version"=>3,"temperature"=>temperature,
        "candidate_filter"=>"parse_then_allowed_APs_then_SAT_before_clustering",
        "allowed_propositions"=>sort!(collect(CONJUNCTION_ALLOWED_PROPOSITIONS)),
        "source_directory"=>abspath(source_directory),
        "manifest_sha256"=>bytes2hex(sha256(JSON3.write(manifest))),"max_k"=>max_k)
    prompt_results = Any[]
    for prompt in manifest["prompts"]
        index = prompt["prompt_index"]
        temperature_label = isnothing(temperature) ? "default" : @sprintf("%.1f",temperature)
        source = joinpath(source_directory,"details",@sprintf("prompt_%03d_T_%s.json",index,temperature_label))
        isfile(source) || error("Missing source samples: $source")
        detail = conjunction_read(source)
        detail["manifest_sha256"] == identity["manifest_sha256"] || error("Source checkpoint manifest mismatch: $source")
        detail["prompt"] == prompt || error("Source checkpoint prompt mismatch: $source")
        detail["temperature"] == temperature || error("Wrong checkpoint temperature")
        samples = detail["samples"]
        [s["repetition"] for s in samples] == collect(1:manifest["repetitions"]) || error("Incomplete or duplicate samples: $source")
        fingerprint = bytes2hex(sha256(read(source)))
        path = joinpath(results_directory,"details",@sprintf("prompt_%03d.json",index))
        if isfile(path)
            result = conjunction_read(path)
            result["analysis_identity"] == identity && result["source_sha256"] == fingerprint || error("Analysis checkpoint changed; use a fresh results_directory")
        else
            excluded = Any[]
            clusters = rank_conjunction_clusters(samples,prompt["ground_truth_ltl"];
                excluded_candidates=excluded)
            result = OrderedDict{String,Any}("analysis_identity"=>identity,
                "source_sha256"=>fingerprint,"prompt_index"=>index,
                "record_id"=>prompt["record_id"],"input_field"=>prompt["input_field"],
                "reference_ltl"=>prompt["ground_truth_ltl"],"samples"=>length(samples),
                "valid_samples"=>count(s -> get(s,"valid",false),samples),
                "satisfiable_samples"=>sum(c["frequency"] for c in clusters;init=0),
                "excluded_candidates"=>excluded,
                "ranked_clusters"=>clusters,
                "conjunctions"=>evaluate_cluster_conjunctions(clusters,prompt["ground_truth_ltl"];max_k=max_k))
            conjunction_save(path,result)
        end
        push!(prompt_results,result)
        println("Conjunction analysis: $(length(prompt_results))/$(length(manifest["prompts"])) prompts; record=$(prompt["record_id"]) clusters=$(length(result["ranked_clusters"]))")
        flush(stdout)
    end
    analysis = merge(identity,OrderedDict{String,Any}("updated_at_utc"=>string(now(UTC)),
        "model"=>get(manifest,"model",nothing),"provider"=>get(manifest,"provider","deepseek"),
        "samples"=>sum(r["samples"] for r in prompt_results),"prompt_results"=>prompt_results))
    report = write_conjunction_tables(results_directory,analysis;
        missing_policy=missing_policy,rescue_definition=rescue_definition)
    println("Table: ",joinpath(results_directory,"comparison.md"))
    return report
end

"""Compute all-prompt SAT-filtered curves from existing outputs, without API calls.

The combined report uses all selected prompts as the denominator. The analysis
subdirectories also retain the existing conditional-on-availability reports.
"""
function analyze_saved_conjunction_curves(; source_directory,
    temperatures=(0.7,1.0), max_k=5,
    results_directory=joinpath(source_directory,"conjunction-curves-sat-filtered"))
    reports = Any[]
    for temperature in temperatures
        label = isnothing(temperature) ? "default" : @sprintf("%.1f",temperature)
        analysis = analyze_cluster_conjunctions(; source_directory=source_directory,
            temperature=temperature, max_k=max_k,
            results_directory=joinpath(results_directory,"analysis-T-$label"))
        report = merge(analysis,OrderedDict{String,Any}(
            "denominator_policy"=>"all_prompts",
            "abstention_definition"=>"No valid, allowed-AP, individually satisfiable candidates; no empty conjunction is evaluated",
            "summary"=>aggregate_all_prompt_conjunctions(analysis["prompt_results"];max_k=max_k)))
        push!(reports,report)
    end
    combined = OrderedDict{String,Any}("schema_version"=>1,
        "updated_at_utc"=>string(now(UTC)),"denominator_policy"=>"all_prompts",
        "missing_policy"=>"saturate","rescue_definition"=>"no_individual_equivalent",
        "temperature_results"=>reports)
    conjunction_save(joinpath(results_directory,"summary.json"),combined)
    metrics = ["logical_good","not_guaranteed","abstain","unsatisfiable","top1_rescue","semantic_rescue"]
    availability = ["available_candidates","available_semantic_clusters","included_representatives","prompts_with_at_least_k_clusters"]
    open(joinpath(results_directory,"summary.csv"),"w") do io
        println(io,join(vcat(["temperature","k","total_prompts","eligible_prompts"],availability,
            metrics,[name*"_count" for name in metrics]),','))
        for report in reports, row in report["summary"]
            println(io,join(vcat(Any[report["temperature"],row["k"],row["total_prompts"],row["eligible_prompts"]],
                Any[row[name] for name in availability],Any[row["rates"][name] for name in metrics],
                Any[row["counts"][name] for name in metrics]),','))
        end
    end
    open(joinpath(results_directory,"comparison.md"),"w") do io
        println(io,"SAT-filtered conjunction curves — $(first(reports)["model"])\n")
        println(io,"Existing outputs only; no model/API calls. Parse, allowed APs, and individual SAT precede semantic clustering. Rank by frequency, breaking ties by earliest sampled member. Conjoin min(k,m) distinct representatives. All probabilities use all selected prompts. With m=0, abstain; no empty conjunction is evaluated.\n")
        for report in reports
            prompts = report["prompt_results"]
            rows = report["summary"]
            excluded = [c for p in prompts for c in p["excluded_candidates"]]
            println(io,"T=$(report["temperature"]): $(length(prompts)) prompts, $(report["samples"]) saved outputs, $(first(rows)["available_candidates"]) retained candidate samples, $(first(rows)["available_semantic_clusters"]) distinct semantic clusters summed across prompts.")
            println(io,"Excluded before clustering: $(count(c -> c["reason"]=="invalid_ltl",excluded)) invalid/incomplete, $(count(c -> c["reason"]=="disallowed_propositions",excluded)) disallowed-AP, $(count(c -> c["reason"]=="unsatisfiable",excluded)) individually UNSAT.\n")
            println(io,"| k | P(good) | P(not guaranteed) | P(abstain) | P(unsat) | P(top-1 rescue) | Semantic rescue | Included representatives (total) | Prompts with ≥k clusters |")
            println(io,"|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for row in rows
                values = [@sprintf("%.1f%%",100*row["rates"][name]) for name in metrics]
                println(io,"| $(row["k"]) | $(join(values," | ")) | $(row["included_representatives"]) | $(row["prompts_with_at_least_k_clusters"]) |")
            end
            println(io)
        end
        println(io,"Semantic rescue: Good(Φ_k) and none of its included representatives is equivalent to the reference. Top-1 rescue: Good(Φ_k) and the top representative is not sufficient. Sufficient means SAT(candidate) and candidate implies reference, over unconstrained infinite traces.\n")
        println(io,"Available candidate samples include duplicate generations after filtering; available semantic clusters are distinct formulas modulo semantic equivalence, counted separately for each prompt. The included-representative count sums min(k,m) over all prompts. Good + not guaranteed + abstain + unsat = 100% in every row.")
    end
    println("All-prompt curves: ",joinpath(results_directory,"comparison.md"))
    return combined
end

if abspath(PROGRAM_FILE)==@__FILE__
    analyze_cluster_conjunctions()
end
