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
    ("./utils/sets/semilinear_set.jl",),
    ("./utils/sets/set_algebra.jl",),
    ("./utils/pclf.jl",),
    ("./symbolic/automaton.jl",),
    ("./system/vector_continuous_system.jl",),
    ("./system/linearization.jl",),
    ("./system/kernel_approximation.jl",),
    ("./system/trajectories.jl",),
    ("./system/transition_synthesis.jl",),
    ("./problem/problems.jl",),
    ("./mapping/list_mapping.jl",),
    ("./mapping/grid.jl",),
    ("./mapping/explicit_grid_mapping.jl",),
    ("./mapping/implicit_grid_mapping.jl",),
    ("./mapping/periodic_grid_mapping.jl",),
    ("./mapping/generic_sets.jl",),
    ("./mapping/abstract_state_set.jl",),
    ("./symbolic/symbolic_model_list.jl",),
    ("./symbolic/allocation.jl",),
    ("./symbolic/multithreading.jl",),
    ("./symbolic/multiprocessing.jl",),
    ("./optim/optimizer_common.jl",),
    ("./optim/LazyEllipsoidAbstraction/lazy_ellipsoid_abstraction.jl", :slow),
    ("./optim/UniformGridAbstraction/unit_test_reachability.jl",),
    ("./optim/UniformGridAbstraction/unit_test_safety.jl",),
    ("./optim/UniformGridAbstraction/unit_test_growth_bound.jl",),
    ("./optim/UniformGridAbstraction/unit_test_linearization.jl",),
    ("./optim/UniformGridAbstraction/uniform_grid_abstraction_reachability.jl",),
    ("./optim/UniformGridAbstraction/uniform_grid_abstraction_safety.jl",),
    ("./optim/UniformGridAbstraction/uniform_grid_abstraction_cosafeLTL.jl",),
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
