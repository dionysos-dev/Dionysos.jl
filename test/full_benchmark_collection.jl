using Test, JSON3, OrderedCollections

# Simulate constants from a pilot already loaded in the user's REPL.
const DEFAULT_INPUT_FIELDS = ["caller_sentinel"]
const PILOT_TEMPERATURES = [0.3, 0.7, 1.0]
include(joinpath(@__DIR__, "..", "FullBenchmarkCollect.jl"))
include(joinpath(@__DIR__, "..", "FullBenchmarkCollect.jl"))
const FB = FullBenchmarkRuntime

function fixture(path)
    records = [Dict{String,Any}(
        "id" => id, "template_id" => id+10, "LTL" => "G(yellow)",
        "task_specific_ap_mapping" => Dict("prop_1" => "yellow"),
    ) for id in 1:2]
    for r in records, field in FB.DEFAULT_INPUT_FIELDS
        # All four fields have identical wording: they must still be separate prompts.
        r[field] = "Always stay in the yellow region. Record $(r["id"])."
    end
    write(path, JSON3.write(records))
end

function response(formula="G(yellow)"; finish_reason="stop")
    return (formula=formula, raw_output="```ltl\n$formula\n```",
        usage=Dict("prompt_tokens" => 100, "completion_tokens" => 7),
        returned_model="deepseek-v4-pro", finish_reason=finish_reason)
end

@testset "Single-file full-benchmark collection, no API calls" begin
    @test DEFAULT_INPUT_FIELDS == ["caller_sentinel"]
    mktempdir() do dir
        dataset = joinpath(dir, "dataset.json")
        output = joinpath(dir, "results.json")
        fixture(dataset)
        settings = (dataset_path=dataset, results_path=output, K=3,
            development_manifest_path=nothing)
        data = prepare_full_benchmark(; settings...)
        @test data["prompt_count"] == 8
        @test data["planned_queries"] == 24
        @test all(isempty(p["samples"]) for p in data["prompts"])
        @test full_benchmark_status(results_path=output)["remaining_queries"] == 24
        @test data["prompts"][1]["translation_prompt"] ==
            FB.build_task_specific_translation_prompt(data["prompts"][1]["source_text"];
                ap_mapping=Dict("prop_1" => "yellow"))

        calls = Any[]
        mock = function(sentence; kwargs...)
            push!(calls, (sentence=sentence, kwargs=kwargs))
            n = length(calls)
            yield() # Exercise interleaved workers without network access.
            formula = ["G(yellow)", "this is invalid LTL", ""][mod1(n, 3)]
            response(formula; finish_reason=isempty(formula) ? "length" : "stop")
        end
        summary = run_full_benchmark(; settings..., request_fn=mock, workers=3,
            checkpoint_every=2, verbose=false)
        @test summary["complete"]
        @test summary["saved_responses"] == length(calls) == 24
        @test summary["total_output_tokens"] == 24*7
        @test summary["mean_output_tokens"] == 7
        @test all(c.kwargs[:temperature] == 0.7 && c.kwargs[:model] == "deepseek-v4-pro" &&
            c.kwargs[:deepseek_thinking_type] == "disabled" && c.kwargs[:return_metadata] &&
            c.kwargs[:ap_mapping] == Dict("prop_1" => "yellow") for c in calls)
        @test all(!haskey(c.kwargs, :ground_truth_ltl) for c in calls)
        @test all(c.sentence in [p["source_text"] for p in data["prompts"]] for c in calls)
        saved = JSON3.read(read(output, String))
        @test all([s["query_index"] for s in p["samples"]] == [1,2,3] for p in saved["prompts"])
        @test any(s["formula"] == "" && s["finish_reason"] == "length" for p in saved["prompts"] for s in p["samples"])
        @test any(s["formula"] == "this is invalid LTL" for p in saved["prompts"] for s in p["samples"])
        @test !isdir(output * ".lock")
        @test Set(readdir(dir)) == Set(["dataset.json", "results.json"])

        # A completed run never requests more responses, even without an API key.
        summary2 = run_full_benchmark(; settings..., request_fn=mock, verbose=false)
        @test length(calls) == 24
        @test summary2 == summary
        @test_throws ErrorException prepare_full_benchmark(; settings..., temperature=1.0)
        @test_throws ErrorException prepare_full_benchmark(; merge(settings, (K=4,))...)
        @test_throws ErrorException prepare_full_benchmark(; settings..., max_tokens=2048)

        # Changing prompt/dataset contents cannot silently reuse old outputs.
        write(dataset, read(dataset, String) * "\n")
        @test_throws ErrorException prepare_full_benchmark(; settings...)
    end

    @testset "Failure persists unsaved responses; manual resume retries only missing query" begin
        mktempdir() do dir
            dataset = joinpath(dir, "dataset.json")
            output = joinpath(dir, "results.json")
            fixture(dataset)
            settings = (dataset_path=dataset, results_path=output, K=3,
                development_manifest_path=nothing)
            attempts = Ref(0)
            fail_once = function(sentence; kwargs...)
                attempts[] += 1
                attempts[] == 5 && error("Mock HTTP 401: invalid API key")
                response()
            end
            @test_throws ErrorException run_full_benchmark(; settings..., request_fn=fail_once,
                workers=1, checkpoint_every=100, verbose=false)
            partial = JSON3.read(read(output, String))
            @test partial["run_state"] == "failed"
            @test partial["summary"]["saved_responses"] == 4
            @test partial["prompts"][2]["errors"][1]["query_index"] == 2
            @test partial["summary"]["failed_attempts"] == 1
            @test !isdir(output * ".lock")
            summary = run_full_benchmark(; settings..., request_fn=fail_once, workers=2, verbose=false)
            @test summary["saved_responses"] == 24
            @test attempts[] == 25 # 24 responses and the one failed request.
            @test summary["failed_attempts"] == 1

            # Corrupted query indices must fail, not cause paid replacement calls.
            damaged = FB.mutable_json(JSON3.read(read(output, String)))
            damaged["prompts"][1]["samples"][2]["query_index"] = 1
            FB.save_results(damaged, output)
            @test_throws ErrorException run_full_benchmark(; settings..., request_fn=fail_once, verbose=false)
            @test attempts[] == 25
            @test !isdir(output * ".lock")
        end
    end

    @testset "Lock prevents a competing writer" begin
        mktempdir() do dir
            output = joinpath(dir, "results.json")
            mkdir(output * ".lock")
            @test_throws ErrorException prepare_full_benchmark(results_path=output)
            @test isdir(output * ".lock") # The rejected process must not remove someone else's lock.
            @test !isfile(output)
        end
    end

    @testset "Save every response; failure retains other in-flight responses" begin
        mktempdir() do dir
            dataset = joinpath(dir, "dataset.json")
            output = joinpath(dir, "results.json")
            fixture(dataset)
            settings = (dataset_path=dataset, results_path=output, K=1,
                development_manifest_path=nothing)
            count = Ref(0)
            inspect_checkpoint = function(sentence; kwargs...)
                # With one worker and checkpoint_every=1, every preceding response is on disk.
                @test full_benchmark_status(results_path=output)["saved_responses"] == count[]
                count[] += 1
                response()
            end
            run_full_benchmark(; settings..., request_fn=inspect_checkpoint,
                workers=1, checkpoint_every=1, verbose=false)
            @test count[] == 8

            output2 = joinpath(dir, "parallel.json")
            count[] = 0
            parallel_failure = function(sentence; kwargs...)
                count[] += 1
                n = count[]
                yield() # Allow all four workers to begin a request before the failure.
                n == 3 && error("Mock credit error")
                yield()
                response()
            end
            @test_throws ErrorException run_full_benchmark(; merge(settings, (results_path=output2,))...,
                request_fn=parallel_failure, workers=4, checkpoint_every=100, verbose=false)
            @test count[] == 4
            status = full_benchmark_status(results_path=output2)
            @test status["saved_responses"] == 3
            @test status["failed_attempts"] == 1
        end
    end

    @testset "Real benchmark plan includes every rewritten field, with development tags" begin
        # Preparing is strictly offline and does not depend on Spot or API credentials.
        mktempdir() do dir
            data = prepare_full_benchmark(results_path=joinpath(dir, "all.json"))
            @test data["record_count"] == 1184
            @test data["prompt_count"] == 4736
            @test data["planned_queries"] == 33152
            @test data["config"]["K"] == 7
            @test data["config"]["atomic_propositions"] == ["blue", "brown", "green", "purple", "yellow"]
            @test length(unique((p["record_id"],p["input_field"]) for p in data["prompts"])) == 4736
            @test all(isempty(p["samples"]) for p in data["prompts"])
            dev_manifest = joinpath(@__DIR__, "..", "Results", "task-specific-temperature-prompt-v2", "manifest.json")
            if isfile(dev_manifest)
                @test count(p -> p["development_prompt_index"] !== nothing, data["prompts"]) == 40
            end
        end
    end
end
