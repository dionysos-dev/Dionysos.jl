using Test, JSON3, SHA, Serialization
include(joinpath(@__DIR__,"worker.jl"))

summary=JSON3.read(read(joinpath(OUTPUT,"summary.json"),String))
@assert summary.complete && length(summary.results)==200
plant=open(deserialize,joinpath(OUTPUT,"plant.jls"))
admissible=Set(admissible_states(plant))
files=unique(String(r.controller_file) for r in summary.results if r.issued)
checks=Any[]
@testset "Exact projected winning regions and controller certificates" begin
    for relative in files
        path=joinpath(OUTPUT,relative)
        c=open(deserialize,path)
        @test !isempty(c.x_win)
        @test Set(c.x_win) ⊆ admissible
        @test Set(c.initial_states)==Set(c.x_win)
        certificate=if c.backend==:buchi
            xs,_=project_buchi_initial_region(plant.ts,plant.labels,c.automaton,c.winning,sort!(collect(admissible)))
            @test xs==c.x_win
            verify_closed_loop_buchi(plant.ts,c.automaton,c.policy,c.next_q,c.x_win)
        else
            @test Set(first.(c.parity_controller.starts))==Set(c.x_win)
            verify_parity_controller(plant.ts,plant.labels,c.parity_controller,EVAL_L2A)
        end
        @test certificate["satisfied"]
        @test c.verification_horizon==:infinite
        @test c.simulation.steps==c.max_simulation_steps==1200
        push!(checks,Dict("controller_file"=>relative,"backend"=>String(c.backend),
            "x_win_size"=>length(c.x_win),"satisfied"=>certificate["satisfied"]))
    end
    for r in summary.results
        @test r.g==Int(r.issued)
        @test r.loss==Int(r.issued && !r.implies_reference)
        if r.region_known
            region=JSON3.read(read(joinpath(OUTPUT,String(r.region_file)),String))
            @test length(region.x_win)==r.x_win_size
            @test r.issued == !isempty(region.x_win)
        else
            @test !r.issued && isnothing(r.x_win_size)
        end
        if r.issued
            @test bytes2hex(open(sha256,joinpath(OUTPUT,String(r.controller_file))))==r.controller_sha256
        end
    end
end
atomic_write_json(joinpath(OUTPUT,"region-controller-audit.json"),Dict(
    "complete"=>true,"distinct_controllers"=>length(files),"admissible_states"=>length(admissible),"results"=>checks))
