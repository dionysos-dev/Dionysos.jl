using Test
using JSON3

# Exercise coexistence with both legacy Main state and the Claude runner.
const PILOT_TEMPERATURES = (0.3, 0.7, 1.0)
const TASK_SPECIFIC_ATOMIC_PROPOSITIONS = ["prop_1"]
request_translation(args...; kwargs...) = error("Stale Main request must not be used")
include(joinpath(@__DIR__, "..", "ClaudeSonnetEvaluate.jl"))
include(joinpath(@__DIR__, "..", "MistralLarge4Evaluate.jl"))

@testset "Mistral loader isolation" begin
    runtime = MistralLarge4Runtime
    include(joinpath(@__DIR__, "..", "MistralLarge4Evaluate.jl"))
    @test runtime === MistralLarge4Runtime
    @test PILOT_TEMPERATURES == (0.3, 0.7, 1.0)
    @test TASK_SPECIFIC_ATOMIC_PROPOSITIONS == ["prop_1"]
    @test runtime.request_translation !== request_translation
    @test runtime.request_translation !== ClaudeSonnetRuntime.request_translation
    @test mistral_large4_settings().temperatures == [0.3, 0.7, 1.0]
    @test sonnet55_settings().temperatures == [nothing]
end

const MISTRAL_TEST_REQUESTS = Any[]
const MISTRAL_TEST_FINISH = Ref("stop")
const MISTRAL_TEST_CONTENT = Ref{Any}("F(blue)")
function MistralLarge4Runtime.post_json(url::String, headers::Vector{Pair{String,String}}, body::String; kwargs...)
    push!(Main.MISTRAL_TEST_REQUESTS, (url=url, headers=Dict(headers), body=JSON3.read(body)))
    parsed = Dict("model"=>MISTRAL_LARGE4_MODEL,
        "choices"=>[Dict("finish_reason"=>Main.MISTRAL_TEST_FINISH[],
            "message"=>Dict("content"=>Main.MISTRAL_TEST_CONTENT[]))],
        "usage"=>Dict("prompt_tokens"=>100, "completion_tokens"=>8, "total_tokens"=>108))
    return (status=200, body=collect(codeunits(JSON3.write(parsed))))
end

@testset "Native Mistral sampling, text chunks, and usage" begin
    ENV["COCO_MISTRAL_TEST_KEY"] = "fixture"
    try
        for temperature in (0.3, 0.7, 1.0)
            response = MistralLarge4Runtime.request_translation("Eventually reach blue.";
                provider=:mistral, model=MISTRAL_LARGE4_MODEL,
                api_url="https://api.mistral.ai/v1/chat/completions",
                api_key_env="COCO_MISTRAL_TEST_KEY", reasoning_effort="none",
                temperature=temperature, return_metadata=true)
            request = last(MISTRAL_TEST_REQUESTS)
            @test request.url == "https://api.mistral.ai/v1/chat/completions"
            @test request.headers["Authorization"] == "Bearer fixture"
            @test !haskey(request.headers,"x-api-key")
            @test request.body["model"] == "mistral-large-4"
            @test request.body["temperature"] == temperature
            @test request.body["reasoning_effort"] == "none"
            @test request.body["max_tokens"] == 1024
            @test !haskey(request.body,"random_seed") && !haskey(request.body,"top_p")
            @test request.body["messages"][1]["content"] == MistralLarge4Runtime.TASK_TRANSLATION_SYSTEM_PROMPT
            @test occursin("AP = {blue, brown, green, purple, yellow}", request.body["messages"][2]["content"])
            @test response.formula == "F(blue)" && response.finish_reason == "stop"
            @test response.returned_model == "mistral-large-4"
            @test MistralLarge4Runtime.score_pilot_sample(response,"F(blue)")["output_tokens"] == 8
        end
        MISTRAL_TEST_CONTENT[] = [Dict("type"=>"thinking", "text"=>"F(yellow)",
            "thinking"=>[Dict("type"=>"text","text"=>"Ignore this reasoning")]),
            Dict("type"=>"text", "text"=>"F(blue)")]
        response = MistralLarge4Runtime.request_translation("Reach blue.";
            provider=:mistral, model=MISTRAL_LARGE4_MODEL, api_key_env="COCO_MISTRAL_TEST_KEY",
            api_url="https://api.mistral.ai/v1/chat/completions", reasoning_effort="none", return_metadata=true)
        @test response.formula == "F(blue)"
        MISTRAL_TEST_FINISH[] = "length"
        truncated = MistralLarge4Runtime.request_translation("Reach blue.";
            provider=:mistral, model=MISTRAL_LARGE4_MODEL, api_key_env="COCO_MISTRAL_TEST_KEY",
            api_url="https://api.mistral.ai/v1/chat/completions", reasoning_effort="none", return_metadata=true)
        @test !MistralLarge4Runtime.score_pilot_sample(truncated,"F(blue)")["valid"]
    finally
        delete!(ENV,"COCO_MISTRAL_TEST_KEY")
        MISTRAL_TEST_CONTENT[] = "F(blue)"
        MISTRAL_TEST_FINISH[] = "stop"
    end
end

@testset "Three-temperature manifest, routing, and resume" begin
    mktempdir() do directory
        runtime = MistralLarge4Runtime
        manifest = prepare_mistral_large4_pilot(; results_directory=directory)
        @test manifest["model"] == "mistral-large-4" && manifest["provider"] == "mistral"
        @test manifest["api_url"] == "https://api.mistral.ai/v1/chat/completions"
        @test manifest["api_key_env"] == "MISTRAL_API_KEY"
        @test manifest["temperatures"] == [0.3, 0.7, 1.0] && manifest["include_temperature"]
        @test manifest["reasoning_effort"] == manifest["thinking_type"] == "none"
        @test manifest["prompt_count"] * manifest["repetitions"] * length(manifest["temperatures"]) == 1200
        @test_throws ArgumentError prepare_mistral_large4_pilot(; results_directory=directory, temperatures=[0.7,0.7])
        @test_throws ArgumentError prepare_mistral_large4_pilot(; results_directory=directory, temperatures=[2.1])
        @test_throws ArgumentError prepare_mistral_large4_pilot(; results_directory=directory, reasoning_effort="low")
        defaults = runtime.prepare_temperature_pilot(; results_directory=joinpath(directory,"native-defaults"),
            provider=:mistral, model=MISTRAL_LARGE4_MODEL, temperatures=[0.3,0.7,1.0])
        @test defaults["api_url"] == manifest["api_url"] && defaults["api_key_env"] == manifest["api_key_env"]
        for prompt in manifest["prompts"], temperature in manifest["temperatures"]
            path = runtime.pilot_detail_path(directory,prompt["prompt_index"],temperature)
            detail = runtime.pilot_detail(manifest,prompt,temperature,path)
            limit = prompt["prompt_index"]==1 ? 9 : 10
            for i in 1:limit
                push!(detail["samples"], Dict("repetition"=>i,"valid"=>true,"satisfiable"=>true,
                    "candidate_formula"=>prompt["ground_truth_ltl"],"equivalent"=>true,
                    "semantically_sufficient"=>true,"output_tokens"=>12,"status"=>"ok"))
            end
            runtime.pilot_save(path,detail)
        end
        called_temperatures = Float64[]
        function fake_request(sentence; kwargs...)
            push!(called_temperatures,kwargs[:temperature])
            @test kwargs[:provider]==:mistral && kwargs[:model]==MISTRAL_LARGE4_MODEL
            @test kwargs[:include_temperature] && kwargs[:reasoning_effort]=="none"
            @test kwargs[:api_url]=="https://api.mistral.ai/v1/chat/completions"
            @test kwargs[:api_key_env]=="MISTRAL_API_KEY"
            prompt = first(manifest["prompts"])
            @test sentence==prompt["source_text"] && kwargs[:ap_mapping]==prompt["task_specific_ap_mapping"]
            return (formula=prompt["ground_truth_ltl"],raw_output=prompt["ground_truth_ltl"],
                usage=Dict("completion_tokens"=>12),returned_model=MISTRAL_LARGE4_MODEL,finish_reason="stop")
        end
        summary = run_mistral_large4_pilot(; results_directory=directory, request_fn=fake_request)
        @test sort(called_temperatures)==[0.3,0.7,1.0]
        @test summary["complete"] && summary["provider"]=="mistral"
        @test length(summary["temperature_results"])==3
        @test all(row["samples"]==400 && row["mean_output_tokens"]==12 && row["per_sample_translation_accuracy"]==1
            for row in summary["temperature_results"])
        table = read(joinpath(directory,"comparison.md"),String)
        @test all(occursin("| $temperature |",table) for temperature in [0.3,0.7,1.0])
        @test occursin("Reasoning effort is fixed at none",table)
        @test run_mistral_large4_pilot(; results_directory=directory,request_fn=fake_request)["complete"]
        @test length(called_temperatures)==3
        @test_throws ErrorException prepare_mistral_large4_pilot(; results_directory=directory,reasoning_effort="high")
    end
end
