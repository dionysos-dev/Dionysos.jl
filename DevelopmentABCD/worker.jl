using JSON3, Serialization, Dates, SHA

const ROOT = abspath(ARGS[1])
const OUTPUT = abspath(ARGS[2])
include(joinpath(ROOT,"evaluate_claude_buchi_synthesis.jl"))
using MathOptInterface, MathematicalSystems
const MOI = MathOptInterface
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const SY = Dionysos.Symbolic
const PB = Dionysos.Problem
const U_BOX = UT.HyperRectangle(SVector(-1.0,-1.0),SVector(1.0,1.0))
const OBSTACLE_MARGIN = 0.12
include(joinpath(ROOT,"System_helper.jl"))

function atomic_serialization(path,value)
    open(path*".tmp","w") do io
        serialize(io,value)
    end
    mv(path*".tmp",path;force=true)
end

function make_plant()
    optimizer = MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer)
    concrete = problem(x1_lb,x1_ub,x2_lb,x2_ub)
    hx = SVector(0.2,0.2,0.3)
    MOI.set(optimizer,MOI.RawOptimizerAttribute("concrete_problem"),concrete)
    MOI.set(optimizer,MOI.RawOptimizerAttribute("state_grid"),DO.GridFree(SVector(0.0,0.0,pi/2),hx))
    MOI.set(optimizer,MOI.RawOptimizerAttribute("input_grid"),DO.GridFree(SVector(0.0,0.0),SVector(0.3,0.3)))
    MOI.set(optimizer,MOI.RawOptimizerAttribute("jacobian_bound"),jacobian_bound())
    MOI.set(optimizer,MOI.RawOptimizerAttribute("time_step"),0.3)
    MOI.set(optimizer,MOI.RawOptimizerAttribute("approx_mode"),Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH)
    MOI.optimize!(optimizer)
    abstract = MOI.get(optimizer,MOI.RawOptimizerAttribute("abstract_system"))
    lo,hi = concrete.region.lb[3],concrete.region.ub[3]
    rect(a,b,c,d)=UT.HyperRectangle(SVector(a,b,lo),SVector(c,d,hi))
    obstacles=[rect(max(0.0,x1_lb[i]-OBSTACLE_MARGIN),max(0.0,x2_lb[i]-OBSTACLE_MARGIN),
        min(10.0,x1_ub[i]+OBSTACLE_MARGIN),min(10.0,x2_ub[i]+OBSTACLE_MARGIN)) for i in eachindex(x1_lb)]
    regions=Dict{Symbol,Any}(:obs=>obstacles,:blue=>rect(6.0,7.5,9.0,9.5),
        :green=>rect(5.0,1.0,7.5,2.5),:purple=>rect(7.0,3.5,8.5,6.5),
        :brown=>rect(4.0,7.0,7.5,8.0),:yellow=>rect(0.5,1.0,2.5,3.5))
    labeler,_=EVAL_L2A.labeler_and_inner_sets(abstract;regions=regions,forbidden_aps=[:obs])
    initial=SVector(4.8,5.0,pi/2)
    box=UT.HyperRectangle(initial-SVector(0.1,0.1,0.1),initial+SVector(0.1,0.1,0.1))
    raw=EVAL_L2A.ts_from_symbolic_onepass(abstract,box;reachable_only=false,x0_abs=nothing)
    cell=UT.HyperRectangle(initial-hx/2,initial+hx/2)
    starts=SY.get_states_from_set(abstract,cell,DO.CENTER)
    isempty(starts) && (starts=SY.get_states_from_set(abstract,cell,DO.OUTER))
    isempty(starts) && error("Cannot map the frozen initial state into the plant")
    ts=EVAL_L2A.TSFromAdj(raw.Adj,[first(starts)])
    labels,_=EVAL_L2A.precompute_labels_and_forbidden(ts,labeler,Symbol[])
    inputs=[collect(SY.get_concrete_input(abstract,i)) for i in 1:EVAL_L2A.ninputs(ts)]
    centers=[collect(SY.get_concrete_state(abstract,i)) for i in 1:EVAL_L2A.nstates(ts)]
    return (;ts,labels,inputs,centers)
end

function logical_sat(formula)
    process=run(pipeline(ignorestatus(`ltlfilt -f $formula --equivalent-to=0 --count`),stdout=devnull))
    process.exitcode in (0,1) || error("Spot SAT check failed")
    return process.exitcode==1
end

function solve_job(job,plant)
    started=time()
    result=Dict{String,Any}("job_id"=>String(job.job_id),"formula"=>String(job.formula),
        "outcome"=>"computational_failure","issued"=>false,"subreason"=>nothing,
        "controller_file"=>nothing,"error"=>nothing)
    stage="logical_satisfiability"
    function phase(name)
        atomic_write_json(joinpath(OUTPUT,"worker-stage.json"),Dict("job_id"=>job.job_id,"stage"=>name))
    end
    try
        phase(stage)
        if !logical_sat(String(job.formula))
            result["outcome"]="logically_inconsistent"
            result["subreason"]="unsatisfiable_conjunction"
            return result
        end
        stage="automaton_translation"; phase(stage)
        aut=EVAL_L2A.translate_ltl_buchi(String(job.formula))
        issubset(Set(aut.aps),TASK_APS_EVAL) || error("Unlabelled propositions in specification")
        result["automaton_states"]=aut.nstates
        # Reading L(x0) is required before testing the initial product state.
        starts=Int[]
        for x0 in EVAL_L2A.initset(plant.ts)
            q=EVAL_L2A.nextstate(aut,aut.init,EVAL_L2A.letter_of(plant.labels[x0],aut.aps))
            if q===nothing
                result["outcome"]="plant_unrealizable"
                result["subreason"]="initial_label_violates_specification"
                return result
            end
            push!(starts,product_id(EVAL_L2A.nstates(plant.ts),x0,q_to_index(aut,q)))
        end
        isempty(starts) && error("Empty configured initial set")
        stage="robust_buchi_game"; phase(stage)
        W,policy,next_q,info=EVAL_L2A.synthesize_buchi_controller(plant.ts,aut,plant.labels;semantics=:robust)
        result["winning_product_states"]=count(W)
        result["total_product_states"]=length(W)
        if !all(W[i] for i in starts)
            result["outcome"]="plant_unrealizable"
            result["subreason"]="initial_product_state_not_winning"
            return result
        end
        stage="universal_closed_loop_verification"; phase(stage)
        certificate=verify_closed_loop_buchi(plant.ts,aut,policy,next_q,EVAL_L2A.initset(plant.ts))
        result["certificate"]=certificate
        certificate["satisfied"] || error("Synthesized policy failed universal infinite-horizon verification")
        stage="controller_export"; phase(stage)
        relative=joinpath("controllers",String(job.job_id)*".jls")
        path=joinpath(OUTPUT,relative)
        mkpath(dirname(path))
        atomic_serialization(path,(;schema_version=1,formula=String(job.formula),
            automaton=aut,policy=policy,next_q=next_q,winning=W,
            initial_product_states=starts,initial_states=EVAL_L2A.initset(plant.ts),
            plant_file="plant.jls",plant_sha256=String(job.plant_sha256),
            semantics=:robust,certificate=certificate))
        result["controller_file"]=relative
        result["controller_sha256"]=bytes2hex(open(sha256,path))
        result["outcome"]="controller_issued"
        result["issued"]=true
        result["subreason"]="verified_and_exported"
    catch err
        err isa InterruptException && rethrow()
        result["error"]=sprint(showerror,err)
        result["subreason"]=occursin("nondeterministic",result["error"]) ? "unsupported_automaton" : stage*"_failure"
    finally
        result["elapsed_seconds"]=time()-started
        result["last_stage"]=stage
    end
    return result
end

function worker_main()
    mkpath(OUTPUT)
    cache=joinpath(OUTPUT,"plant.jls")
    plant=if isfile(cache)
        open(deserialize,cache)
    else
        value=make_plant()
        atomic_serialization(cache,value)
        value
    end
    # Compile the game and verifier before the per-formula wall-clock budget.
    tiny=EVAL_L2A.TSFromAdj([[[1]]],[1])
    aut=EVAL_L2A.Automaton(["blue"],0,1,true,Dict{Tuple{Int,UInt64},Int}((0,UInt64(1))=>0),BitSet([0]))
    _,p,n,_=EVAL_L2A.synthesize_buchi_controller(tiny,aut,[Set([:blue])];semantics=:robust)
    verify_closed_loop_buchi(tiny,aut,p,n,[1])["satisfied"] || error("Verifier warm-up failed")
    atomic_write_json(joinpath(OUTPUT,"worker-ready.json"),Dict("pid"=>getpid(),
        "plant_states"=>EVAL_L2A.nstates(plant.ts),"plant_inputs"=>EVAL_L2A.ninputs(plant.ts),
        "initial_states"=>EVAL_L2A.initset(plant.ts),
        "initial_labels"=>[sort!(string.(collect(plant.labels[i]))) for i in EVAL_L2A.initset(plant.ts)]))
    for line in eachline(stdin)
        isempty(strip(line)) && continue
        request=JSON3.read(line)
        job=JSON3.read(read(String(request.job_file),String))
        result=solve_job(job,plant)
        atomic_write_json(String(request.result_file),result)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    worker_main()
end
