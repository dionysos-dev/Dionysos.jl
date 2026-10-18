# ==============================================================================
# Control synthesis and animation of the normal-form biped experiment, from the saved abstraction
# ==============================================================================
#
# Runs `examples/BipedRobot/biped_4d_normal_form.jl` on the abstraction file written by
# `merge_chunks.jl` (`ABSTRACTION_FILE`): nothing is computed but the synthesis, the closed loop and the
# figures, which are written to `DIONYSOS_OUTPUT_DIR`. The file is never recomputed here: if it is
# missing, this stops (it would be a sequential computation of the whole abstraction).
#
# Rerun it as often as needed with other starts or tolerances: they do not enter the abstraction.
# Put them in a Julia file and point `DIONYSOS_SOLVE_PARAMS` at it; it is included before the example,
# so it sets the variables the example reads, e.g.
#
#     FOOTHOLD_TOLERANCE = 0.03
#     INITIAL_VELOCITY = (-0.44, -0.36)
#     Q0 = SVector(0.0, 0.5, -0.4, -0.9)      # StaticArrays is loaded for it
#     TARGET_VELOCITY = (-0.3, -0.1)
#
# (`FOOTHOLD_TOLERANCE`, `Q0`, `INITIAL_VELOCITY`, `TARGET_VELOCITY`.) What defines the abstraction
# (`DISC`, `X_BOUNDS`, `ONE_JOINT_PER_STEP`, `OBSTACLE`) is the example's: the file is refused if it differs.
#
# Environment:
#   DIONYSOS_ABSTRACTION_FILE   the saved abstraction (default: slurm_array/out/normal_form_abstraction.jld2)
#   DIONYSOS_SOLVE_PARAMS       optional Julia file with the variables above
#   DIONYSOS_SCENARIO           name of the outputs (default: job<SLURM_JOB_ID>)
#   DIONYSOS_OUTPUT_DIR         directory of the outputs (default: slurm_array/out/results)

import Pkg
Pkg.instantiate()

using StaticArrays   # for the params file

const EXAMPLE = abspath(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "..",
        "examples",
        "BipedRobot",
        "biped_4d_normal_form.jl",
    ),
)

abstraction_path = abspath(
    get(
        ENV,
        "DIONYSOS_ABSTRACTION_FILE",
        joinpath(@__DIR__, "out", "normal_form_abstraction.jld2"),
    ),
)
isfile(abstraction_path) || error(
    "$abstraction_path does not exist: run the chunks and the merge first (submit_all.sh), " *
    "or set DIONYSOS_ABSTRACTION_FILE",
)

# Read by the example before it defaults them.
ABSTRACTION_FILE = abstraction_path
SCENARIO_NAME =
    Symbol(get(ENV, "DIONYSOS_SCENARIO", "job" * get(ENV, "SLURM_JOB_ID", "local")))
OUTPUT_DIR = abspath(get(ENV, "DIONYSOS_OUTPUT_DIR", joinpath(@__DIR__, "out", "results")))
RUN_PLOTS = true

params = get(ENV, "DIONYSOS_SOLVE_PARAMS", "")
if !isempty(params)
    isfile(params) || error("DIONYSOS_SOLVE_PARAMS = $params does not exist")
    include(abspath(params))      # may override the variables above too (e.g. SCENARIO_NAME)
end

mkpath(OUTPUT_DIR)
prefix = joinpath(OUTPUT_DIR, "biped_4d_normal_form_$(SCENARIO_NAME)")
isempty(params) || cp(params, prefix * "_params.jl"; force = true)   # a record of what produced it

@info "Control synthesis" abstraction_path SCENARIO_NAME OUTPUT_DIR params
t = @elapsed include(EXAMPLE)

if isfile(prefix * ".gif")
    @info "Done" seconds = round(t; digits = 1) gif = prefix * ".gif"
else
    @error "No animation was written: the controller has no action at the start state, or the " *
           "closed loop did not run (see the lines above). Try a larger FOOTHOLD_TOLERANCE or " *
           "another start." prefix
    exit(2)
end
