using Test, JSON3, SHA
include(joinpath(@__DIR__, "..", "FullBenchmarkExtend.jl"))
include(joinpath(@__DIR__, "..", "FullBenchmarkExtend.jl"))
const EFB=FullBenchmarkRuntime

function extension_fixture(path)
    records=[Dict{String,Any}("id"=>i,"template_id"=>i,"LTL"=>"G(yellow)",
        "task_specific_ap_mapping"=>Dict("prop_1"=>"yellow")) for i in 1:2]
    for r in records, field in EFB.DEFAULT_INPUT_FIELDS
        r[field]="Stay in yellow. Record $(r["id"])."
    end
    write(path,JSON3.write(records))
end

extension_response(text)=(formula=text,raw_output=text,usage=Dict("completion_tokens"=>4),
    returned_model="fixture",finish_reason="stop")

@testset "Extend a saved run without repeating or replacing old responses" begin
    mktempdir() do dir
        dataset=joinpath(dir,"dataset.json")
        source=joinpath(dir,"R2.json")
        target=joinpath(dir,"R4.json")
        extension_fixture(dataset)
        old_calls=Ref(0)
        old_request=function(sentence;kwargs...)
            old_calls[]+=1
            extension_response("old response $(old_calls[])")
        end
        settings=(dataset_path=dataset,development_manifest_path=nothing)
        EFB.run_full_benchmark(;settings...,results_path=source,K=2,
            request_fn=old_request,workers=2,verbose=false)
        old=EFB.mutable_json(JSON3.read(read(source,String)))
        push!(old["prompts"][1]["errors"],Dict("query_index"=>1,"error"=>"Old failed attempt","at_utc"=>"fixture"))
        EFB.checkpoint!(old,source)
        original_hash=sha256(read(source))
        @test old_calls[]==16
        extended=extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)
        @test extended["config"]["K"]==extended["config"]["repetitions"]==4
        @test extended["planned_queries"]==32
        @test extended["summary"]["saved_responses"]==16
        @test extended["summary"]["remaining_queries"]==16
        @test extended["summary"]["completed_prompts"]==0
        @test extended["summary"]["failed_attempts"]==1
        @test all(p["samples"]==o["samples"] && p["errors"]==o["errors"] for (p,o) in zip(extended["prompts"],old["prompts"]))
        @test extended["extension_history"][1]["source_sha256"]==bytes2hex(original_hash)
        @test sha256(read(source))==original_hash

        # An API failure during the extension must still resume at the first missing new index.
        new_calls=Ref(0)
        new_request=function(sentence;kwargs...)
            new_calls[]+=1
            new_calls[]==3 && error("Mock balance error")
            extension_response("new response $(new_calls[])")
        end
        @test_throws ErrorException EFB.run_full_benchmark(;settings...,results_path=target,K=4,
            request_fn=new_request,workers=1,checkpoint_every=50,verbose=false)
        partial=extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)
        @test partial["summary"]["saved_responses"]==18
        @test partial["prompts"][1]["samples"][3]["query_index"]==3
        summary=EFB.run_full_benchmark(;settings...,results_path=target,K=4,
            request_fn=new_request,workers=2,verbose=false)
        @test summary["complete"] && summary["saved_responses"]==32
        @test new_calls[]==17 # 16 additional responses plus one failed request.
        done=EFB.mutable_json(JSON3.read(read(target,String)))
        @test all(p["samples"][1:2]==o["samples"] for (p,o) in zip(done["prompts"],old["prompts"]))
        @test all([s["query_index"] for s in p["samples"]]==[1,2,3,4] for p in done["prompts"])
        @test sha256(read(source))==original_hash
        @test extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)["summary"]["remaining_queries"]==0
        @test !isdir(source*".lock") && !isdir(target*".lock")
        @test_throws ArgumentError extend_full_benchmark(;settings...,source_path=source,results_path=source,K=4)
        @test_throws ArgumentError extend_full_benchmark(;settings...,source_path=source,results_path=target,K=2)

        # Reject a destination with altered historical samples without overwriting it.
        done["prompts"][1]["samples"][1]["raw_output"]="changed historical output"
        EFB.save_results(done,target)
        damaged_hash=sha256(read(target))
        @test_throws ErrorException extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)
        @test sha256(read(target))==damaged_hash
        @test sha256(read(source))==original_hash
    end

    @testset "Incomplete source and changed experiment are rejected" begin
        mktempdir() do dir
            dataset=joinpath(dir,"dataset.json")
            source=joinpath(dir,"R2.json")
            target=joinpath(dir,"R4.json")
            extension_fixture(dataset)
            settings=(dataset_path=dataset,development_manifest_path=nothing)
            EFB.prepare_full_benchmark(;settings...,results_path=source,K=2)
            @test_throws ErrorException extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)
            @test !isfile(target)
            @test !isdir(source*".lock") && !isdir(target*".lock")
            # A changed benchmark must not silently extend an old experiment.
            write(dataset,read(dataset,String)*"\n")
            @test_throws ErrorException extend_full_benchmark(;settings...,source_path=source,results_path=target,K=4)
            @test !isfile(target)
        end
    end
end
