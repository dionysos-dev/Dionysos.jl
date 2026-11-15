using Test
include(joinpath(@__DIR__,"worker.jl"))

@testset "Issuance requires an exported verified controller" begin
    ts=EVAL_L2A.TSFromAdj([[[1]]],[1])
    plant=(;ts,labels=[Set([:blue])])
    job(formula,id)=(;formula,job_id=id,plant_sha256="fixture")
    issued=solve_job(job("G(blue)","issued"),plant)
    @test issued["outcome"]=="controller_issued" && issued["issued"]
    artifact=open(deserialize,joinpath(OUTPUT,issued["controller_file"]))
    @test artifact.certificate["satisfied"]
    @test artifact.policy==[1]
    @test all(artifact.winning)
    @test bytes2hex(open(sha256,joinpath(OUTPUT,issued["controller_file"])))==issued["controller_sha256"]
    impossible=solve_job(job("blue & !blue","unsat"),plant)
    @test impossible["outcome"]=="logically_inconsistent" && !impossible["issued"]
    losing=solve_job(job("G(yellow)","losing"),plant)
    @test losing["outcome"]=="plant_unrealizable" && !losing["issued"]
    unsupported=solve_job(job("F(G(blue))","unsupported"),plant)
    @test unsupported["outcome"]=="computational_failure" && !unsupported["issued"]
    @test unsupported["subreason"]=="unsupported_automaton"
end

@testset "An accepting prefix cannot certify infinite acceptance" begin
    ts=EVAL_L2A.TSFromAdj([[[2]],[[2]]],[1])
    labels=[Set([:blue]),Set{Symbol}()]
    aut=EVAL_L2A.Automaton(["blue"],0,2,true,
        Dict{Tuple{Int,UInt64},Int}((0,UInt64(1))=>0,(0,UInt64(0))=>1,
            (1,UInt64(1))=>0,(1,UInt64(0))=>1),BitSet([0]))
    next_q=EVAL_L2A._precompute_next_q(aut,labels)
    certificate=verify_closed_loop_buchi(ts,aut,[1,1,1,1],next_q,[1])
    @test !certificate["satisfied"]
    @test certificate["has_reachable_nonaccepting_cycle"]
end
