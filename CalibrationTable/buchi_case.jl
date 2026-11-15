# Reuse the production fixed-point solver and independent infinite-run verifier.
# Parity is dispatched separately by run.py, with its explicit larger limits.
include(joinpath(ARGS[1], "InitialRegionABCD", "worker.jl"))
job = JSON3.read(read(ARGS[4], String))
plant = open(deserialize, String(job.plant_cache))
admissible = admissible_states(plant)
length(admissible) == 44643 || error("Admissible plant domain changed")
# Fail loudly if dispatch disagrees; never silently use the old parity caps.
translate_deterministic_buchi(String(job.formula), EVAL_L2A)
result = solve_region_job(job, plant, String[] , admissible)
result["backend"] == "parity" && error("Unexpected parity fallback in Buchi-only dispatch")
atomic_write_json(joinpath(OUTPUT, "verified-result.json"), result)
