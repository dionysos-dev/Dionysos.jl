using Test
include(joinpath(@__DIR__,"..","ClusterConjunction.jl"))

function conjunction_fixture(formulas)
    return [Dict("candidate_formula"=>f,"valid"=>!isempty(f),"repetition"=>i)
        for (i,f) in enumerate(formulas)]
end

@testset "Native Claude default-temperature checkpoints" begin
    mktempdir() do source
        prompts = [OrderedDict{String,Any}("prompt_index"=>i,"record_id"=>i,
            "input_field"=>"fixture","ground_truth_ltl"=>"F(blue)") for i in 1:2]
        manifest = OrderedDict{String,Any}("prompt_version"=>2,
            "temperatures"=>[nothing],"include_temperature"=>false,
            "atomic_propositions"=>sort!(collect(CONJUNCTION_ALLOWED_PROPOSITIONS)),
            "model"=>"claude-sonnet-5-5","provider"=>"anthropic",
            "repetitions"=>2,"prompts"=>prompts)
        conjunction_save(joinpath(source,"manifest.json"),manifest)
        manifest_hash = bytes2hex(sha256(JSON3.write(manifest)))
        for prompt in prompts
            formulas = prompt["prompt_index"]==1 ? ["G(blue)","F(prop_1)"] : ["",""]
            conjunction_save(joinpath(source,"details",@sprintf("prompt_%03d_T_default.json",prompt["prompt_index"])),
                OrderedDict{String,Any}("manifest_sha256"=>manifest_hash,"prompt"=>prompt,
                    "temperature"=>nothing,"samples"=>conjunction_fixture(formulas)))
        end
        report = analyze_cluster_conjunctions(;source_directory=source,temperature=nothing)
        @test report["temperature"]===nothing && report["model"]=="claude-sonnet-5-5"
        @test report["samples"]==4
        @test all(r["eligible_prompts"]==1 && r["total_prompts"]==2 for r in report["summary"])
        @test all(r["rates"]["logical_good"]==1 && r["rates"]["semantic_rescue"]==1 for r in report["summary"])
        table = joinpath(source,"conjunction-temperature-default-validated","comparison.md")
        @test occursin("claude-sonnet-5-5, provider-default temperature",read(table,String))
        @test analyze_cluster_conjunctions(;source_directory=source,temperature=nothing)["summary"]==report["summary"]
        @test_throws ErrorException analyze_cluster_conjunctions(;source_directory=source,temperature=0.7)
        curves = analyze_saved_conjunction_curves(;source_directory=source,
            temperatures=(nothing,),max_k=5,results_directory=joinpath(source,"all-prompt-curves"))
        @test curves["denominator_policy"]=="all_prompts"
        @test all(r["rates"]["logical_good"]==0.5 && r["rates"]["abstain"]==0.5
            for r in only(curves["temperature_results"])["summary"])
        @test isfile(joinpath(source,"all-prompt-curves","summary.csv"))
    end
end

@testset "All-prompt curves preserve abstention and candidate counts" begin
    pools = [["G(blue)"],["blue","blue","yellow"],["purple","!purple"],String[]]
    references = ["F(blue)","blue & yellow","green","blue"]
    prompts = Any[]
    for (pool,reference) in zip(pools,references)
        clusters = rank_conjunction_clusters(conjunction_fixture(pool),reference)
        push!(prompts,Dict("ranked_clusters"=>clusters,
            "satisfiable_samples"=>sum(c["frequency"] for c in clusters;init=0),
            "conjunctions"=>evaluate_cluster_conjunctions(clusters,reference;max_k=5)))
    end
    rows = aggregate_all_prompt_conjunctions(prompts)
    @test all(r["total_prompts"]==4 && r["eligible_prompts"]==3 for r in rows)
    @test all(r["rates"]["abstain"]==0.25 for r in rows)
    @test all(r["available_candidates"]==6 && r["available_semantic_clusters"]==5 for r in rows)
    @test rows[1]["rates"]["logical_good"]==0.25
    @test rows[1]["rates"]["not_guaranteed"]==0.5
    @test rows[1]["rates"]["unsatisfiable"]==0
    @test rows[2]["rates"]["logical_good"]==0.5
    @test rows[2]["rates"]["not_guaranteed"]==0
    @test rows[2]["rates"]["unsatisfiable"]==0.25
    @test rows[2]["rates"]["top1_rescue"]==0.25
    @test rows[2]["rates"]["semantic_rescue"]==0.5
    @test [r["included_representatives"] for r in rows]==[3,5,5,5,5]
    @test [r["prompts_with_at_least_k_clusters"] for r in rows]==[3,2,0,0,0]
    @test all(sum(r["rates"][name] for name in ("logical_good","not_guaranteed","abstain","unsatisfiable"))==1 for r in rows)
    empty = aggregate_all_prompt_conjunctions([last(prompts)])
    @test all(r["rates"]["abstain"]==1 && r["rates"]["logical_good"]==0 for r in empty)
end

@testset "Top-1 versus synergistic rescue" begin
    # The first two formulas are insufficient alone, but sufficient together.
    # Rank 3 is sufficient by itself, so including it ends synergistic rescue.
    ranked=rank_conjunction_clusters(conjunction_fixture(
        ["blue","blue","blue","yellow","yellow","blue & yellow"]),"blue & yellow")
    prefixes=evaluate_cluster_conjunctions(ranked,"blue & yellow")
    rows=aggregate_cluster_conjunctions([Dict("conjunctions"=>prefixes)])
    @test rows[1]["counts"]["top1_rescue"]==rows[1]["counts"]["synergistic_rescue"]==0
    @test rows[2]["counts"]["top1_rescue"]==rows[2]["counts"]["synergistic_rescue"]==1
    @test rows[3]["counts"]["top1_rescue"]==1
    @test rows[3]["counts"]["synergistic_rescue"]==0
    @test rows[3]["counts"]["logical_good"]==1
    @test rows[3]["rates"]==rows[10]["rates"]
    # Selection rescue can occur by adding an already-sufficient lower cluster.
    selected=rank_conjunction_clusters(conjunction_fixture(
        ["F(yellow)","F(yellow)","G(blue)"]),"F(blue)")
    metrics=aggregate_cluster_conjunctions([Dict("conjunctions"=>evaluate_cluster_conjunctions(selected,"F(blue)"))])
    @test metrics[2]["counts"]["top1_rescue"]==1
    @test metrics[2]["counts"]["synergistic_rescue"]==0
    @test all(r["counts"]["synergistic_rescue"]<=r["counts"]["top1_rescue"]<=r["counts"]["logical_good"] for r in metrics)
end

@testset "Explicit parse and AP checks precede SAT and equivalence" begin
    excluded=Any[]
    samples=conjunction_fixture(["F(","F(prop_1)","F(obstacle)","blue & !blue","F(blue)"])
    ranked=rank_conjunction_clusters(samples,"F(blue)";excluded_candidates=excluded)
    @test length(ranked)==1 && first(ranked)["sample_indices"]==[5]
    @test count(c -> c["reason"]=="invalid_ltl",excluded)==1
    @test count(c -> c["reason"]=="disallowed_propositions",excluded)==2
    @test count(c -> c["reason"]=="unsatisfiable",excluded)==1
    @test !haskey(CONJUNCTION_SAT_CACHE,"F(prop_1)")
    @test !haskey(CONJUNCTION_SAT_CACHE,"F(obstacle)")
    @test all(!("F(prop_1)" in key) && !("F(obstacle)" in key) for key in keys(CONJUNCTION_EQUIVALENCE_CACHE))
    @test length(filter_conjunction_candidates([Dict("candidate_formula"=>"G(blue)","repetition"=>1)]))==1
    @test length(filter_conjunction_candidates(conjunction_fixture(["G(\"blue\")"])))==1
    @test length(filter_conjunction_candidates(conjunction_fixture(["true"])))==1
end

@testset "Frequency ranking and equivalence" begin
    ranked=rank_conjunction_clusters(conjunction_fixture(
        ["F(blue)","G(yellow)","!(G(!blue))","G(yellow)",""]),"F(blue)")
    @test length(ranked)==2
    @test [c["frequency"] for c in ranked]==[2,2]
    @test first(ranked)["sample_indices"]==[1,3]
    @test first(ranked)["equivalent_to_reference"]
    @test !last(ranked)["equivalent_to_reference"]
end

@testset "Conjunction rescue and vacuous implication" begin
    ranked=rank_conjunction_clusters(conjunction_fixture(["blue","blue","yellow"]),"blue & yellow")
    rows=evaluate_cluster_conjunctions(ranked,"blue & yellow")
    @test rows[1]["not_guaranteed"] && !rows[1]["logical_good"]
    @test rows[2]["satisfiable"] && rows[2]["logical_good"]
    @test rows[2]["semantic_rescue_no_individual_sufficient"]
    @test rows[2]["semantic_rescue_no_individual_equivalent"]
    @test conjunction_equivalent(rows[2]["phi_k"],"blue & yellow")
    impossible=rank_conjunction_clusters(conjunction_fixture(["blue","blue","!blue"]),"yellow")
    unsat=evaluate_cluster_conjunctions(impossible,"yellow")
    @test unsat[2]["unsatisfiable"] && unsat[2]["implies_reference"]
    @test !unsat[2]["logical_good"] && !unsat[2]["not_guaranteed"]
    @test !unsat[2]["semantic_rescue_no_individual_equivalent"]
    stronger=rank_conjunction_clusters(conjunction_fixture(["G(blue)"]),"F(blue)")
    strong=evaluate_cluster_conjunctions(stronger,"F(blue)")
    @test strong[1]["semantic_rescue_no_individual_equivalent"]
    @test !strong[1]["semantic_rescue_no_individual_sufficient"]
    contradiction=rank_conjunction_clusters(conjunction_fixture(["blue & !blue"]),"yellow")
    @test isempty(contradiction)
    @test isempty(evaluate_cluster_conjunctions(contradiction,"yellow"))
    @test_throws ErrorException conjunction_spot_match("F(","--equivalent-to=0")
end

@testset "Candidate SAT filtering happens before clustering and ranking" begin
    excluded=Any[]
    samples=conjunction_fixture(["blue & !blue","blue & !blue","blue",
        "yellow & !yellow","blue & !blue",""])
    ranked=rank_conjunction_clusters(samples,"blue";excluded_candidates=excluded)
    @test length(ranked)==1
    @test first(ranked)["representative"]=="blue"
    @test first(ranked)["frequency"]==1
    @test first(ranked)["sample_indices"]==[3]
    @test all(c["satisfiable"] for c in ranked)
    @test count(c -> c["reason"]=="unsatisfiable",excluded)==4
    @test count(c -> c["reason"]=="invalid_ltl",excluded)==1
    @test !any(any(occursin("blue & !blue",f) || occursin("yellow & !yellow",f) for f in key)
        for key in keys(CONJUNCTION_EQUIVALENCE_CACHE))
end

@testset "Saturation, denominators, and exhaustive categories" begin
    a=rank_conjunction_clusters(conjunction_fixture(["G(blue)"]),"F(blue)")
    b=rank_conjunction_clusters(conjunction_fixture(["blue","!blue"]),"yellow")
    prompts=[Dict("conjunctions"=>evaluate_cluster_conjunctions(a,"F(blue)")),
        Dict("conjunctions"=>evaluate_cluster_conjunctions(b,"yellow"))]
    rows=aggregate_cluster_conjunctions(prompts)
    @test length(rows)==10
    @test all(r["eligible_prompts"]==2 for r in rows)
    @test rows[1]["rates"]["satisfiable"]==1
    @test rows[2]["rates"]["unsatisfiable"]==0.5
    @test rows[2]["rates"]["implies_reference"]==1
    @test rows[2]["rates"]["logical_good"]==0.5
    @test rows[2]["rates"]["semantic_rescue"]==0.5
    @test rows[2]["rates"]==rows[10]["rates"]
    @test all(sum(r["rates"][k] for k in ["logical_good","not_guaranteed","unsatisfiable"])==1 for r in rows)
    strict=aggregate_cluster_conjunctions(prompts;missing_policy="strict")
    @test strict[1]["eligible_prompts"]==2 && strict[2]["eligible_prompts"]==1
    @test strict[3]["eligible_prompts"]==0 && strict[3]["rates"]["satisfiable"]===nothing
    empty=aggregate_cluster_conjunctions([Dict("conjunctions"=>Any[])])
    @test empty[1]["eligible_prompts"]==0
    @test_throws ArgumentError aggregate_cluster_conjunctions(prompts;missing_policy="invent")
end
