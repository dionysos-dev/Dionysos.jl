# Test driver.
#
# Fast dev loop: `Pkg.test(; test_args = ["--fast"])` skips suites tagged `:slow` — slow end-to-end
# optimization-solver tests whose coverage isn't worth the wall-clock while iterating. A plain
# `Pkg.test()` and CI run everything. Each file's run time is printed, with a slowest-first summary at
# the end, so slow suites are easy to spot and tag `:slow`.
const FAST_TESTS = "--fast" in ARGS

# (path, tags...). Tag a suite `:slow` to exclude it from `--fast`.
const TEST_FILES = [
    ("./aqua.jl", :slow),  # quality/meta gate (ambiguities, stale deps): run pre-commit / in CI
    ("./utils/data_structures/sorted_vector_set.jl",),
    ("./utils/data_structures/tree.jl",),
    ("./utils/functions.jl",),
    ("./utils/numeric/scalar_optimization.jl",),
    ("./utils/sets/rectangle.jl",),
    ("./utils/sets/ellipsoid.jl",),
    ("./utils/sets/ellipsoid_intersection.jl",),
    ("./utils/sets/semilinear_set.jl",),
    ("./utils/sets/set_algebra.jl",),
    ("./utils/periodic.jl",),
    ("./utils/pclf.jl",),
    ("./symbolic/automaton.jl",),
    ("./system/vector_continuous_system.jl",),
    ("./system/affine_approximation.jl",),
    ("./system/approximation.jl",),
    ("./system/controllers.jl",),
    ("./system/pid_controller.jl",),
    ("./system/controller_serialization.jl",),
    ("./system/trajectories.jl",),
    ("./system/transition_synthesis.jl",),
    ("./problem/problems.jl",),
    ("./mapping/list_mapping.jl",),
    ("./mapping/grid.jl",),
    ("./mapping/explicit_grid_mapping.jl",),
    ("./mapping/implicit_grid_mapping.jl",),
    ("./mapping/periodic_grid_mapping.jl",),
    ("./mapping/multi_level_mapping.jl",),
    ("./mapping/generic_sets.jl",),
    ("./mapping/abstract_state_set.jl",),
    ("./symbolic/symbolic_model_list.jl",),
    ("./symbolic/clock_lift.jl",),
    ("./symbolic/lift_per_slice.jl",),
    ("./symbolic/states_satisfying.jl",),
    ("./symbolic/allocation.jl",),
    ("./symbolic/multithreading.jl",),
    ("./symbolic/multiprocessing.jl",),
    ("./optim/optimizer_common.jl",),
    ("./optim/jump_frontend.jl",),  # canonical JuMP entry (∂/final on Dionysos.Optimizer)
    ("./optim/LazyEllipsoidAbstraction/lazy_ellipsoid_abstraction.jl", :slow),
    # Direct discrete-automaton controller synthesis (no abstraction build).
    ("./optim/discrete_systems/reachability.jl",),
    ("./optim/discrete_systems/safety.jl",),
    # UniformGridAbstraction: abstraction-build modes + end-to-end control specs.
    ("./optim/UniformGridAbstraction/growth_bound.jl",),
    ("./optim/UniformGridAbstraction/linearized.jl",),
    ("./optim/UniformGridAbstraction/user_defined.jl",),
    ("./optim/UniformGridAbstraction/reachability.jl",),
    ("./optim/UniformGridAbstraction/safety.jl",),
    ("./optim/UniformGridAbstraction/reach_and_stay.jl",),
    ("./optim/UniformGridAbstraction/cosafe_ltl.jl",),
    ("./optim/UniformGridAbstraction/clock_lifted_continuous.jl",),
    # Trajectory generators / certifiers and the PCLF bisimulation quotient.
    ("./optim/trajectory_generators.jl",),
    ("./optim/trajectory_certifiers.jl", :slow),
    ("./optim/PCLFBisimulationQuotient/bisimulation_quotient.jl", :slow),
    ("./optim/UniformEllipsoidAbstraction/uniform_ellipsoid_abstraction.jl", :slow),
    ("./optim/HybridSystemAbstraction/hybrid_system_abstraction.jl", :slow),
    ("./optim/GolLazarBelta/gol_lazar_belta.jl", :slow),
    ("./optim/NonLinear/non_linear.jl", :slow),
    ("./regression/uniform_grid_abstraction.jl", :slow),  # golden-output regression net
]

const _timings = Tuple{String, Float64}[]
for entry in TEST_FILES
    path = entry[1]
    tags = entry[2:end]
    if FAST_TESTS && (:slow in tags)
        println("[skip] ", path)
        continue
    end
    dt = @elapsed include(path)
    println("[time] ", lpad(string(round(dt; digits = 1)), 7), "s  ", path)
    push!(_timings, (path, dt))
end

println("\n===== Test file timings (slowest first) =====")
for (path, sec) in sort(_timings; by = x -> -x[2])
    println(lpad(string(round(sec; digits = 1)), 8), "s  ", path)
end
println(
    "total: ",
    round(sum(last, _timings; init = 0.0); digits = 1),
    "s over ",
    length(_timings),
    " files",
)
