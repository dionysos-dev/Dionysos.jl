using Test, JSON3
include(joinpath(@__DIR__,"..","VLTLEvaluate.jl"))
include(joinpath(@__DIR__,"..","VLTLEvaluate.jl"))
const V=VLTLRuntime

response(f;finish_reason="stop")=(formula=f,raw_output=f,usage=Dict("prompt_tokens"=>10,"completion_tokens"=>3),
    returned_model="deepseek-v4-pro",finish_reason=finish_reason)

@testset "VLTL prompt, collection, resuming, and semantic analysis" begin
    benchmark=V.vltl_plan()
    @test benchmark["prompt_count"]==42
    @test benchmark["planned_queries"]==420
    @test benchmark["config"]["atomic_propositions"]==["prop_1","prop_2","prop_3","prop_4"]
    @test all(occursin("AP = {prop_1, prop_2, prop_3, prop_4}",p["translation_prompt"]) for p in benchmark["prompts"])
    @test all(!occursin(r"\b(blue|brown|green|purple|yellow)\b",p["translation_prompt"]) for p in benchmark["prompts"])
    payload=V.request_body(benchmark["prompts"][1]["translation_prompt"],benchmark["config"])
    @test payload["messages"][2]["content"]==benchmark["prompts"][1]["translation_prompt"]
    @test payload["model"]=="deepseek-v4-pro" && payload["temperature"]==0.7 && payload["thinking"]["type"]=="disabled"
    @test !haskey(payload,"ground_truth_ltl")
    # Generic demonstrations are fixed; benchmark references never enter the request.
    mktempdir() do dir
        dataset=joinpath(dir,"fixture.json");path=joinpath(dir,"raw.json")
        records=[Dict("NL"=>"prop_1 eventually holds; prop_2 may hold.","LTL"=>"F(prop_1)"),
                 Dict("NL"=>"Both prop_1 and prop_2 hold now.","LTL"=>"prop_1 & prop_2"),
                 Dict("NL"=>"prop_1 holds now.","LTL"=>"prop_1"),
                 Dict("NL"=>"Eventually prop_2 holds.","LTL"=>"F(prop_2)")]
        write(dataset,JSON3.write(records))
        settings=(dataset_path=dataset,results_path=path,R=3)
        prepared=prepare_vltl(;settings...)
        @test isempty(prepared["prompts"][1]["samples"])
        calls=Ref(0)
        pools=[["G(prop_1)","G(prop_1)","G(prop_1)"],
               ["prop_1","prop_1","prop_2"],
               ["prop_1","!prop_1","prop_1"],
               ["blue","prop_2 & !prop_2",""]]
        mock=function(prompt,config)
            calls[]+=1;index=findfirst(p->p["translation_prompt"]==prompt,prepared["prompts"])
            j=mod1(calls[],3)
            @test config["temperature"]==0.7 && config["thinking_type"]=="disabled"
            @test config["atomic_propositions"]==["prop_1","prop_2"]
            return response(pools[index][j])
        end
        summary=run_vltl(;settings...,workers=1,request_fn=mock,verbose=false)
        @test summary["complete"] && calls[]==12
        @test summary["mean_output_tokens"]==3
        @test run_vltl(;settings...,request_fn=(args...)->error("Must reuse completed collection"),verbose=false)==summary
        @test !isdir(path*".lock")
        @test_throws ErrorException prepare_vltl(;settings...,temperature=1.0)
        analysis=analyze_vltl(source_path=path,max_k=3)
        @test analysis["complete"]
        @test [r["rates"]["logical_good"] for r in analysis["summary"]]==[0.5,0.5,0.5]
        @test [r["rates"]["unsatisfiable"] for r in analysis["summary"]]==[0,0.25,0.25]
        @test all(r["rates"]["abstain"]==0.25 for r in analysis["summary"])
        @test analysis["exclusion_counts"]==Dict("invalid_ltl"=>1,"disallowed_propositions"=>1,"unsatisfiable"=>1)
        resumed=analyze_vltl(source_path=path,max_k=3)
        @test resumed["summary"]==analysis["summary"]
        @test !isdir(splitext(path)[1]*"-conjunctions/summary.json.lock")
    end
    mktempdir() do dir
        dataset=joinpath(dir,"fixture.json");path=joinpath(dir,"raw.json")
        write(dataset,JSON3.write([Dict("NL"=>"Eventually prop_1 holds.","LTL"=>"F prop_1")]))
        settings=(dataset_path=dataset,results_path=path,R=3)
        calls=Ref(0)
        fail=function(args...)
            calls[]+=1
            calls[]==2 && error("fixture API failure")
            response("F prop_1")
        end
        @test_throws ErrorException run_vltl(;settings...,workers=1,request_fn=fail,verbose=false)
        @test vltl_status(results_path=path)["saved_responses"]==1
        @test !isdir(path*".lock")
        remaining=Ref(0)
        success=(args...)->(remaining[]+=1;response("F prop_1"))
        @test run_vltl(;settings...,workers=2,request_fn=success,verbose=false)["complete"]
        @test remaining[]==2
    end
end
