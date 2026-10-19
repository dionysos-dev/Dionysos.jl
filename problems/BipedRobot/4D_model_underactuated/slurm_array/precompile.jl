# ==============================================================================
# Precompile the packages of the normal-form biped experiment, on a compute node
# ==============================================================================
#
# The tasks of the array load Dionysos, Symbolics, Plots, ... each in its own Julia process. If the
# precompilation cache of the node they run on is not valid (never built, or built by another kind of
# CPU: the cache holds machine code and is rejected on another CPU model), every task rebuilds it before
# doing any work, minutes of CPU each, with all the tasks fighting for the same lock files. This job builds
# it once, before the array, and the array waits for it (`submit_all.sh`).
#
# It compiles only the packages the experiment loads (not the whole project: MeshCat, RigidBodyDynamics,
# ... are never used by it), then checks, in a fresh Julia process, that loading them is fast. The check
# is the useful output: if it says the load is slow, the cache is rejected on this node, see README §6.
#
# Environment: JULIA_CPU_TARGET (if set, the same value must be exported for all the jobs: sbatch passes
# the environment of `submit_all.sh` on), DIONYSOS_LOAD_LIMIT (seconds, default 120).

ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"     # `instantiate` must not precompile the whole project
import Pkg
using Dates

stamp(msg) = @info "[precompile $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"

# What the experiment loads (examples/BipedRobot/biped_4d_normal_form.jl and the files it includes).
const PACKAGES = [
    "StaticArrays",
    "MathematicalSystems",
    "LazySets",
    "MathOptInterface",
    "Dionysos",
    "Symbolics",
    "IntervalArithmetic",
    "JLD2",
    "Plots",
]

stamp(
    "start: host $(gethostname()), CPU $(Sys.CPU_NAME) ($(Sys.CPU_THREADS) threads), " *
    "Julia $(VERSION), JULIA_CPU_TARGET=$(get(ENV, "JULIA_CPU_TARGET", "(default)")), " *
    "depot $(first(DEPOT_PATH)), precompile tasks $(get(ENV, "JULIA_NUM_PRECOMPILE_TASKS", "(default)"))",
)
t = @elapsed Pkg.instantiate()
stamp("instantiate: $(round(t; digits = 1)) s")
t = @elapsed Pkg.precompile(PACKAGES)
stamp("precompile of the packages of the experiment: $(round(t; digits = 1)) s")

# Loading them is what the tasks do. The first load compiles what Pkg.precompile leaves out (the package
# extensions that need several packages together, e.g. Dionysos with Symbolics and IntervalArithmetic);
# the second one is the measure: it must be fast, since nothing is left to compile.
load = `$(Base.julia_cmd()) --project=$(Base.active_project()) -e "using $(join(PACKAGES, ", "))"`
t1 = @elapsed run(load)
stamp("first load of the packages in a fresh process: $(round(t1; digits = 1)) s")
t2 = @elapsed run(load)
stamp("second load: $(round(t2; digits = 1)) s")

limit = parse(Float64, get(ENV, "DIONYSOS_LOAD_LIMIT", "120"))
if t2 > limit
    stamp(
        "PRECOMPILE FAILED: loading the packages still takes $(round(t2; digits = 1)) s (more than " *
        "$limit s) after precompiling: the cache is rejected on this node, so the tasks would " *
        "recompile it too. Another CPU than the one that built it, or another JULIA_CPU_TARGET / " *
        "Julia: see README §6.",
    )
    exit(1)
end
stamp("PRECOMPILE OK: loading takes $(round(t2; digits = 1)) s")
