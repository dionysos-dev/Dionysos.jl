using Test
include(joinpath(@__DIR__,"worker.jl"))

function fixture(formula,id,ts,labels)
    plant=(;ts,labels)
    admissible=admissible_states(plant)
    candidates=EVAL_L2A.TSFromAdj(ts.Adj,admissible)
    aps=export_parity_plant(candidates,labels,PARITY_PLANT,EVAL_L2A)
    return solve_region_job((;formula,job_id=id,plant_sha256="fixture"),plant,aps,admissible;max_steps=37)
end

@testset "Initial-region projection and existential issuance" begin
    ts=EVAL_L2A.TSFromAdj([[[1]],[[2]]],[2])
    labels=[Set([:blue]),Set{Symbol}()]
    result=fixture("G blue","buchi-region",ts,labels)
    @test result["issued"] && result["backend"]=="buchi"
    @test result["x_win_size"]==1
    artifact=open(deserialize,joinpath(OUTPUT,result["controller_file"]))
    @test artifact.x_win==[1] && artifact.initial_states==[1]
    @test artifact.admissible_initial_states==[1,2]
    @test artifact.simulation.steps==37
    @test artifact.certificate["satisfied"]
    result=fixture("F G blue","parity-region",ts,labels)
    @test result["issued"] && result["backend"]=="parity"
    @test result["x_win_size"]==1
    artifact=open(deserialize,joinpath(OUTPUT,result["controller_file"]))
    @test artifact.x_win==[1]
    @test artifact.parity_controller.starts[1][1]==1
    @test artifact.certificate["satisfied"]
    # The accepting automaton sink wins as a raw product state, but cannot be
    # initialized from an all-not-blue plant. It must not produce issuance.
    blank=EVAL_L2A.TSFromAdj([[[1]]],[1])
    result=fixture("F blue","raw-product-is-insufficient",blank,[Set{Symbol}()])
    @test result["winning_product_states"]>0
    @test result["x_win_size"]==0 && !result["issued"]
    @test result["outcome"]=="plant_unrealizable"
    # Next is evaluated after the current label, not at the current state.
    next=EVAL_L2A.TSFromAdj([[[2]],[[2]]],[1])
    result=fixture("X blue","next",next,labels)
    @test result["x_win_size"]==0 && !result["issued"]
    result=fixture("blue & !blue","unsat",ts,labels)
    @test result["region_known"] && result["x_win_size"]==0
    @test result["outcome"]=="logically_inconsistent"
end
