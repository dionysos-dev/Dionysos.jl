using Test
if !isdefined(@__MODULE__, :FullBenchmarkConjunctionRuntime)
    include(joinpath(@__DIR__,"..","FullBenchmarkConjunction.jl"))
end
const FBC=FullBenchmarkConjunctionRuntime

function full_fixture(path)
    pools=[["G(blue)","G(blue)","G(blue)"], ["blue","blue","yellow"],
        ["purple","!purple","purple"], ["F(prop_1)","","blue & !blue"]]
    references=["F(blue)","blue & yellow","green","F(blue)"]
    prompts=Any[]
    for (i,(pool,reference)) in enumerate(zip(pools,references))
        push!(prompts,Dict("prompt_index"=>i,"record_id"=>i,"input_field"=>"fixture",
            "source_text"=>"fixture $i","ground_truth_ltl"=>reference,
            "samples"=>[Dict("query_index"=>j,"formula"=>f,"status"=>"returned","finish_reason"=>"stop")
                for (j,f) in enumerate(pool)]))
    end
    data=Dict("stage"=>"LLM_collection_only","config"=>Dict("prompt_version"=>2,
        "atomic_propositions"=>["blue","brown","green","purple","yellow"],
        "K"=>3,"model"=>"fixture","temperature"=>0.7),
        "prompts"=>prompts,"prompt_count"=>4,"planned_queries"=>12)
    write(path,FBC.JSON3.write(data))
end

@testset "Full-file conjunction adapter matches development definitions" begin
    mktempdir() do dir
        source=joinpath(dir,"source.json")
        full_fixture(source)
        output=joinpath(dir,"analysis")
        report=analyze_full_benchmark_conjunctions(source_path=source,
            results_directory=output,max_k=3,workers=3,checkpoint_every=1)
        rows=report["summary"]
        @test report["complete"] && report["completed_prompts"]==4
        @test [r["rates"]["logical_good"] for r in rows]==[0.25,0.5,0.5]
        @test [r["rates"]["not_guaranteed"] for r in rows]==[0.5,0,0]
        @test all(r["rates"]["abstain"]==0.25 for r in rows)
        @test [r["rates"]["unsatisfiable"] for r in rows]==[0,0.25,0.25]
        @test [r["rates"]["top1_rescue"] for r in rows]==[0,0.25,0.25]
        @test [r["rates"]["semantic_rescue"] for r in rows]==[0.25,0.5,0.5]
        @test [r["included_representatives"] for r in rows]==[3,5,5]
        @test all(r["available_candidates"]==9 && r["available_semantic_clusters"]==5 for r in rows)
        @test report["exclusion_counts"]==Dict("invalid_ltl"=>1,"disallowed_propositions"=>1,"unsatisfiable"=>1)
        @test isfile(joinpath(output,"summary.json")) && isfile(joinpath(output,"comparison.md"))
        @test !isdir(joinpath(output,".analysis.lock"))
        @test FBC.JSON3.read(read(source,String))["planned_queries"]==12

        # Resuming skips all completed logical analyses, even with empty Spot caches.
        empty!(FBC.CONJUNCTION_SAT_CACHE)
        empty!(FBC.CONJUNCTION_IMPLIES_CACHE)
        empty!(FBC.CONJUNCTION_EQUIVALENCE_CACHE)
        resumed=analyze_full_benchmark_conjunctions(source_path=source,
            results_directory=output,max_k=3,workers=1)
        @test resumed["summary"]==rows
        @test isempty(FBC.CONJUNCTION_SAT_CACHE) && isempty(FBC.CONJUNCTION_IMPLIES_CACHE) && isempty(FBC.CONJUNCTION_EQUIVALENCE_CACHE)
        @test_throws ErrorException analyze_full_benchmark_conjunctions(source_path=source,
            results_directory=output,max_k=4)
        @test !isdir(joinpath(output,".analysis.lock"))
        # Any source change invalidates a resumed analysis, never silently merges it.
        write(source,read(source,String)*"\n")
        @test_throws ErrorException analyze_full_benchmark_conjunctions(source_path=source,
            results_directory=output,max_k=3)
    end
    mktempdir() do dir
        source=joinpath(dir,"source.json")
        full_fixture(source)
        data=FBC.conjunction_read(source)
        pop!(data["prompts"][1]["samples"])
        FBC.conjunction_save(source,data)
        @test_throws ErrorException analyze_full_benchmark_conjunctions(source_path=source,
            results_directory=joinpath(dir,"analysis"))
        @test !isdir(joinpath(dir,"analysis"))
    end
end
