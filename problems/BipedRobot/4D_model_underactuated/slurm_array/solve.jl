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

# `instantiate` must not precompile the whole project (MeshCat, RigidBodyDynamics, ... that the experiment
# never loads): the precompilation is the job of `precompile.jl`, once, before the array.
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
import Pkg
Pkg.instantiate()

using StaticArrays   # for the params file
using Dates
import JLD2

stamp(msg) = @info "[solve $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
peak_gb() = round(Sys.maxrss() / 2^30; digits = 1)

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

try
    stamp("start: job $(get(ENV, "SLURM_JOB_ID", "-")), memory requested " *
          "$(get(ENV, "SLURM_MEM_PER_NODE", "?")) MB, Julia $(VERSION), $(gethostname())")
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
    n = JLD2.jldopen(f -> haskey(f, "n_transitions") ? f["n_transitions"] : nothing, abstraction_path, "r")
    stamp("abstraction $abstraction_path: $(round(filesize(abstraction_path) / 2^30; digits = 1)) GB, " *
          "$(n === nothing ? "?" : n) transitions (in memory, an automaton takes of the order of " *
          "$(n === nothing ? "?" : round(Int, 200 * n / 2^30)) GB for them)")

    # Read by the example before it defaults them.
    global ABSTRACTION_FILE = abstraction_path
    global SCENARIO_NAME =
        Symbol(get(ENV, "DIONYSOS_SCENARIO", "job" * get(ENV, "SLURM_JOB_ID", "local")))
    global OUTPUT_DIR =
        abspath(get(ENV, "DIONYSOS_OUTPUT_DIR", joinpath(@__DIR__, "out", "results")))
    global RUN_PLOTS = true

    params = get(ENV, "DIONYSOS_SOLVE_PARAMS", "")
    if !isempty(params)
        isfile(params) || error("DIONYSOS_SOLVE_PARAMS = $params does not exist")
        include(abspath(params))      # may override the variables above too (e.g. SCENARIO_NAME)
    end

    mkpath(OUTPUT_DIR)
    prefix = joinpath(OUTPUT_DIR, "biped_4d_normal_form_$(SCENARIO_NAME)")
    isempty(params) || cp(params, prefix * "_params.jl"; force = true)   # a record of what produced it

    stamp("setting up the model, then loading the abstraction; scenario $SCENARIO_NAME, " *
          "outputs $prefix*, params '$params'")
    t = @elapsed include(EXAMPLE)
    stamp("synthesis and closed loop done in $(round(t; digits = 1)) s (peak memory $(peak_gb()) GB)")

    if isfile(prefix * ".gif")
        stamp("SOLVE OK: $(prefix * ".gif")")
    else
        stamp("SOLVE FAILED: no animation was written: the controller has no action at the start " *
              "state, or the closed loop did not run (see the lines above). Try a larger " *
              "FOOTHOLD_TOLERANCE or another start.")
        exit(2)
    end
catch e
    stamp("SOLVE FAILED after peak memory $(peak_gb()) GB: " * (e isa OutOfMemoryError ?
        "OutOfMemoryError: raise --mem of submit_solve.sbatch" : sprint(showerror, e)))
    @error "SOLVE FAILED" exception = (e, catch_backtrace())
    exit(1)
end
