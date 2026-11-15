# Versioned separately from the completed fixed-start experiment.
include(joinpath(@__DIR__, "..", "DevelopmentABCD", "worker.jl"))
include(joinpath(@__DIR__, "ParityBackend.jl"))
using .ParityBackend

const MAX_STEPS = length(ARGS)>=3 ? parse(Int,ARGS[3]) : 1200
const BRIDGE = joinpath(@__DIR__, "build", "spot_parity_bridge")
const PARITY_PLANT = joinpath(OUTPUT, "plant.parity.bin")

function admissible_states(plant)
    # This abstraction already represents X_free. Retain an explicit obstacle
    # check so a differently labelled cache cannot silently admit obstacles.
    return [x for x in 1:EVAL_L2A.nstates(plant.ts) if !(:obs in plant.labels[x])]
end

function project_buchi_initial_region(ts, labels, aut, winning, admissible)
    xs, products = Int[], Int[]
    for x in admissible
        q = EVAL_L2A.nextstate(aut, aut.init, EVAL_L2A.letter_of(labels[x], aut.aps))
        isnothing(q) && continue
        pid = product_id(EVAL_L2A.nstates(ts), x, q_to_index(aut,q))
        if winning[pid]
            push!(xs,x);push!(products,pid)
        end
    end
    return xs,products
end

function solve_region_job(job, plant, aps, admissible; max_steps=MAX_STEPS)
    started=time()
    result=Dict{String,Any}("job_id"=>String(job.job_id), "formula"=>String(job.formula),
        "outcome"=>"computational_failure", "issued"=>false, "subreason"=>nothing,
        "controller_file"=>nothing, "error"=>nothing, "backend"=>nothing,
        "region_known"=>false, "x_win_size"=>nothing, "region_file"=>nothing,
        "admissible_initial_state_count"=>length(admissible),
        "max_simulation_steps"=>max_steps, "finite_trace_is_proof"=>false)
    stage="logical_satisfiability"
    phase(name)=atomic_write_json(joinpath(OUTPUT,"worker-stage.json"),Dict("job_id"=>job.job_id,"stage"=>name))
    function record_region(xs)
        xs=sort!(unique(xs))
        relative=joinpath("regions",String(job.job_id)*".json")
        mkpath(joinpath(OUTPUT,"regions"))
        atomic_write_json(joinpath(OUTPUT,relative),Dict("job_id"=>job.job_id,
            "formula"=>job.formula,"x_win"=>xs,"size"=>length(xs),
            "admissible_initial_state_count"=>length(admissible),
            "projection"=>"(x, delta(q_init, label(x))) belongs to the product winning set"))
        result["region_known"]=true;result["x_win_size"]=length(xs)
        result["region_file"]=relative
        result["region_sha256"]=bytes2hex(open(sha256,joinpath(OUTPUT,relative)))
        return xs
    end
    try
        phase(stage)
        if !logical_sat(String(job.formula))
            record_region(Int[])
            result["outcome"]="logically_inconsistent";result["subreason"]="unsatisfiable_conjunction"
            return result
        end
        stage="buchi_translation";phase(stage)
        aut=try
            translate_deterministic_buchi(String(job.formula),EVAL_L2A)
        catch err
            err isa BuchiUnavailable || rethrow()
            result["buchi_unavailable_reason"]=sprint(showerror,err)
            nothing
        end
        payload=if !isnothing(aut)
            result["backend"]="buchi"
            issubset(Set(aut.aps),TASK_APS_EVAL) || error("Unlabelled specification proposition")
            result["automaton_states"]=aut.nstates
            stage="global_buchi_game";phase(stage)
            W,policy,next_q,_=EVAL_L2A.synthesize_buchi_controller(plant.ts,aut,plant.labels;semantics=:robust)
            result["winning_product_states"]=count(W)
            result["total_product_states"]=length(W)
            xs,products=project_buchi_initial_region(plant.ts,plant.labels,aut,W,admissible)
            record_region(xs)
            if isempty(xs)
                result["outcome"]="plant_unrealizable";result["subreason"]="empty_admissible_initial_winning_region"
                return result
            end
            stage="universal_buchi_verification";phase(stage)
            certificate=verify_closed_loop_buchi(plant.ts,aut,policy,next_q,xs)
            certificate["satisfied"] || error("Büchi verification failed for the projected initial region")
            stage="bounded_simulation";phase(stage)
            states,memory,inputs,status=EVAL_L2A.simulate_buchi_policy(plant.ts,aut,plant.labels,policy,next_q;
                x0=first(xs),q0=aut.init,max_steps=max_steps,stop_at_accept=false,consume_initial_label=true)
            status==:max_steps || error("Büchi simulation failed: $status")
            simulation=(;states,memory,inputs,status="horizon_reached",steps=length(inputs))
            (;backend=:buchi,automaton=aut,policy,next_q,winning=W,
                initial_product_states=products,initial_states=xs,certificate,simulation)
        else
            result["backend"]="parity"
            stage="global_spot_parity_game";phase(stage)
            parity=synthesize_parity(String(job.formula),plant.ts,plant.labels,EVAL_L2A;
                bridge=BRIDGE,plant_binary=PARITY_PLANT,aps=aps,
                directory=joinpath(OUTPUT,"parity-jobs",String(job.job_id)),max_steps=max_steps)
            result["automaton_states"]=parity.info.automaton_states
            result["winning_product_states"]=parity.info.winning_product_states
            result["game_nodes"]=parity.info.game_nodes;result["game_edges"]=parity.info.game_edges
            result["parity_priorities"]=parity.info.priorities
            xs=record_region(Int.(parity.info.winning_initial_states))
            @assert parity.info.winning == !isempty(xs)
            if isempty(xs)
                result["outcome"]="plant_unrealizable";result["subreason"]="empty_admissible_initial_winning_region"
                return result
            end
            @assert Set(first.(parity.controller.starts))==Set(xs)
            stage="universal_parity_verification";phase(stage)
            certificate=verify_parity_controller(plant.ts,plant.labels,parity.controller,EVAL_L2A)
            certificate["satisfied"] || error("Parity verification failed for the projected initial region")
            stage="bounded_simulation";phase(stage)
            simulation=simulate_parity_controller(plant.ts,plant.labels,parity.controller,EVAL_L2A;max_steps=max_steps)
            (;backend=:parity,parity_controller=parity.controller,initial_states=xs,
                certificate,simulation,specification_hoa=read(parity.hoa,String))
        end
        stage="controller_export";phase(stage)
        relative=joinpath("controllers",String(job.job_id)*".jls");mkpath(joinpath(OUTPUT,"controllers"))
        artifact=merge(payload,(;schema_version=3,formula=String(job.formula),
            plant_file="plant.jls",plant_sha256=String(job.plant_sha256),semantics=:robust,
            admissible_initial_states=admissible,x_win=payload.initial_states,
            max_simulation_steps=max_steps,verification_horizon=:infinite))
        atomic_serialization(joinpath(OUTPUT,relative),artifact)
        result["controller_file"]=relative
        result["controller_sha256"]=bytes2hex(open(sha256,joinpath(OUTPUT,relative)))
        result["certificate"]=payload.certificate
        result["simulation_steps"]=payload.simulation.steps;result["simulation_status"]="horizon_reached"
        result["outcome"]="controller_issued";result["issued"]=true;result["subreason"]="verified_nonempty_initial_region"
    catch err
        err isa InterruptException && rethrow()
        result["error"]=sprint(showerror,err)
        result["subreason"]=occursin("resource limit",result["error"]) ? "resource_limit" : stage*"_failure"
    finally
        result["elapsed_seconds"]=time()-started;result["last_stage"]=stage
    end
    return result
end

function region_main()
    MAX_STEPS>0 || error("Invalid simulation horizon")
    mkpath(OUTPUT)
    cache=joinpath(OUTPUT,"plant.jls")
    plant=if isfile(cache)
        open(deserialize,cache)
    else
        p=make_plant();atomic_serialization(cache,p);p
    end
    admissible=admissible_states(plant)
    candidate_ts=EVAL_L2A.TSFromAdj(plant.ts.Adj,admissible)
    aps=export_parity_plant(candidate_ts,plant.labels,PARITY_PLANT,EVAL_L2A)
    tiny=EVAL_L2A.TSFromAdj([[[1]]],[1])
    aut=EVAL_L2A.Automaton(["blue"],0,1,true,Dict{Tuple{Int,UInt64},Int}((0,UInt64(1))=>0),BitSet([0]))
    _,p,n,_=EVAL_L2A.synthesize_buchi_controller(tiny,aut,[Set([:blue])];semantics=:robust)
    verify_closed_loop_buchi(tiny,aut,p,n,[1])["satisfied"] || error("Warm-up failed")
    atomic_write_json(joinpath(OUTPUT,"admissible-initial-states.json"),Dict("states"=>admissible,"size"=>length(admissible)))
    atomic_write_json(joinpath(OUTPUT,"worker-ready.json"),Dict("pid"=>getpid(),
        "plant_states"=>EVAL_L2A.nstates(plant.ts),"plant_inputs"=>EVAL_L2A.ninputs(plant.ts),
        "admissible_initial_state_count"=>length(admissible),"initial_domain"=>"all_admissible_plant_states",
        "max_simulation_steps"=>MAX_STEPS,"cached_plant_initial_states_ignored"=>EVAL_L2A.initset(plant.ts)))
    for line in eachline(stdin)
        isempty(strip(line)) && continue
        request=JSON3.read(line);job=JSON3.read(read(String(request.job_file),String))
        atomic_write_json(String(request.result_file),solve_region_job(job,plant,aps,admissible))
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    region_main()
end
