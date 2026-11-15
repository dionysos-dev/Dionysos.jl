using JSON3, Serialization, SHA, Test
const ROOT=abspath(ARGS[1])
const OUTPUT=abspath(ARGS[2])
include(joinpath(ROOT,"evaluate_claude_buchi_synthesis.jl"))

summary=JSON3.read(read(joinpath(OUTPUT,"summary.json"),String))
@assert summary.complete
@assert length(summary.results)==200
plant=open(deserialize,joinpath(OUTPUT,"plant.jls"))
plant_hash=bytes2hex(open(sha256,joinpath(OUTPUT,"plant.jls")))
issued=filter(r -> r.issued,summary.results)
files=unique(String(r.controller_file) for r in issued)
@testset "Independent verification of every issued controller artifact" begin
    for file in files
        controller=open(deserialize,joinpath(OUTPUT,file))
        @test controller.plant_sha256==plant_hash
        @test controller.semantics==:robust
        @test all(controller.winning[i] for i in controller.initial_product_states)
        certificate=verify_closed_loop_buchi(plant.ts,controller.automaton,
            controller.policy,controller.next_q,controller.initial_states)
        @test certificate["satisfied"]
    end
    for row in summary.results
        @test row.issued==(row.outcome=="controller_issued")
        if row.issued
            @test bytes2hex(open(sha256,joinpath(OUTPUT,String(row.controller_file))))==row.controller_sha256
            controller=open(deserialize,joinpath(OUTPUT,String(row.controller_file)))
            @test controller.formula==row.formula
        end
    end
end
println("Independently re-verified ",length(files)," controller artifacts for ",length(issued)," issued prompt/budget cases.")
