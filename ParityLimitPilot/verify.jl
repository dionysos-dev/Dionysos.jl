# Independent production verification of a native bridge strategy.
include(joinpath(ARGS[1],"InitialRegionABCD","worker.jl"))
plan=JSON3.read(read(ARGS[4],String))
info=JSON3.read(read(joinpath(OUTPUT,"spot-result.json"),String))
plant=open(deserialize,String(plan.plant_cache))
aps=open(String(plan.plant_binary),"r") do io
    String(read(io,8))=="ABCDTS01" || error("Bad native plant")
    nx,nu,ni,na=[Int(ltoh(read(io,UInt32))) for _ in 1:4]
    nx==EVAL_L2A.nstates(plant.ts) || error("Plant dimension mismatch")
    nu==EVAL_L2A.ninputs(plant.ts) || error("Plant input mismatch")
    [String(read(io,Int(ltoh(read(io,UInt32))))) for _ in 1:na]
end
controller=ParityBackend.read_policy(joinpath(OUTPUT,"spot-policy.bin"),aps)
xs=sort!(unique(Int.(info.winning_initial_states)))
info.winning && info.verified && !isempty(xs) || error("Native strategy was not verified winning")
Set(first.(controller.starts))==Set(xs) || error("Native initial region mismatch")
admissible=admissible_states(plant)
issubset(Set(xs),Set(admissible)) || error("Winning initial region contains inadmissible states")
certificate=verify_parity_controller(plant.ts,plant.labels,controller,EVAL_L2A)
certificate["satisfied"] || error("Independent Julia parity verification failed")
simulation=simulate_parity_controller(plant.ts,plant.labels,controller,EVAL_L2A;max_steps=MAX_STEPS)
artifact=(;schema_version=3,formula=String(plan.formula),backend=:parity,
    parity_controller=controller,initial_states=xs,x_win=xs,admissible_initial_states=admissible,
    certificate,simulation,specification_hoa=read(String(plan.hoa_path),String),
    plant_file=String(plan.plant_cache),plant_sha256=String(plan.plant_sha256),
    semantics=:robust,max_simulation_steps=MAX_STEPS,verification_horizon=:infinite)
atomic_serialization(joinpath(OUTPUT,"controller.jls"),artifact)
atomic_write_json(joinpath(OUTPUT,"region.json"),Dict("formula"=>plan.formula,"x_win"=>xs,"size"=>length(xs)))
result=Dict("issued"=>true,"outcome"=>"controller_issued","certificate"=>certificate,
    "x_win_size"=>length(xs),"region_known"=>true,"simulation_steps"=>simulation.steps,
    "simulation_status"=>simulation.status,"verification_horizon"=>"infinite",
    "controller_file"=>joinpath(OUTPUT,"controller.jls"),
    "controller_sha256"=>bytes2hex(open(sha256,joinpath(OUTPUT,"controller.jls"))),
    "region_file"=>joinpath(OUTPUT,"region.json"),
    "region_sha256"=>bytes2hex(open(sha256,joinpath(OUTPUT,"region.json"))))
atomic_write_json(joinpath(OUTPUT,"verified-result.json"),result)
