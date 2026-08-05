# Cross-solver benchmark cases (current API).
#
# Each case is a NamedTuple `(; method, problem, run)` where `run()` performs the full
# setup-and-solve for one (solver family, problem) pair and returns the configured optimizer.
# The runner in `solver_suite.jl` times `run()` end-to-end (setup + abstraction build + synthesis)
# with a per-case try/catch, so one failing case never aborts the suite.
#
# Setups mirror the corresponding tests under `test/optim/**`, so they track the live API.

import LazySets
import Random
using StaticArrays
using LinearAlgebra
using JuMP
import MathOptInterface as MOI

using Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const OP = Dionysos.Optim
const AB = OP.Abstraction

# Open optimizers only (no license-gated solvers): OSQP/HiGHS/Ipopt/Pavito for the MIQP
# families, Clarabel for the SDP (ellipsoid) families, CDDLib for polyhedra.
import CDDLib
import Clarabel
import HiGHS
import Ipopt
import OSQP
import Pavito

# ---------------------------------------------------------------------------
# Problem library (subfoldered under ../problems); each include defines a module.
# ---------------------------------------------------------------------------
const PROBLEMS_DIR = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS_DIR, "DCDC", "dcdc_converter.jl"))
include(joinpath(PROBLEMS_DIR, "PathPlanning", "path_planning.jl"))
include(joinpath(PROBLEMS_DIR, "GolLazarBelta", "gol_lazar_belta.jl"))
include(joinpath(PROBLEMS_DIR, "PwaSystem", "pwa_system.jl"))
include(joinpath(PROBLEMS_DIR, "NonLinear", "non_linear.jl"))

# ---------------------------------------------------------------------------
# Shared solver instances (silent).
# ---------------------------------------------------------------------------
const MIP_SOLVER = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
const QP_SOLVER = optimizer_with_attributes(
    OSQP.Optimizer,
    "polish" => 1,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true,
)
const CONT_SOLVER = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
const MIQP_SOLVER = optimizer_with_attributes(
    Pavito.Optimizer,
    "mip_solver" => MIP_SOLVER,
    "cont_solver" => CONT_SOLVER,
    MOI.Silent() => true,
)
const SDP_SOLVER = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# ---------------------------------------------------------------------------
# Case runners
# ---------------------------------------------------------------------------

# UniformGridAbstraction on the DC-DC converter (safety, GROWTH mode).
function run_uniform_grid_dcdc()
    concrete_problem = DCDC.problem()
    state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(2.0 / 4.0e3, 2.0 / 4.0e3))
    input_grid = MP.GridFree(SVector(1), SVector(1))

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
    MOI.optimize!(optimizer)
    return optimizer
end

# UniformGridAbstraction on path planning (reach-avoid, GROWTH mode).
function run_uniform_grid_path_planning()
    concrete_problem = PathPlanning.problem(; simple = true)
    state_grid = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2))
    input_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3))

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        PathPlanning.jacobian_bound(),
    )
    MOI.optimize!(optimizer)
    return optimizer
end

# A feasible depth-9 Gol–Lazar–Belta instance shared by the MIQP families.
_glb_problem() = GolLazarBelta.problem(
    CDDLib.Library(),
    Float64;
    N = 9,
    q_0 = 8,
    x_0 = [1.5, -2.5],
    zero_cost = false,
)

# Bemporad–Morari (MIQP) on Gol–Lazar–Belta.
function run_bemporad_morari_glb()
    optimizer = MOI.instantiate(
        optimizer_with_attributes(
            OP.BemporadMorari.Optimizer{Float64},
            "continuous_solver" => QP_SOLVER,
            "mixed_integer_solver" => MIQP_SOLVER,
            "indicator" => false,
            "print_level" => 0,
        ),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), _glb_problem())
    MOI.optimize!(optimizer)
    return optimizer
end

# Branch and bound (with the Hybrid Dual Dynamic Programming lower bound) on Gol–Lazar–Belta.
function run_branch_and_bound_glb()
    optimizer = MOI.instantiate(
        optimizer_with_attributes(
            OP.BranchAndBound.Optimizer{Float64},
            "continuous_solver" => QP_SOLVER,
            "mixed_integer_solver" => MIQP_SOLVER,
            "max_iter" => 1000,
            "max_time" => 300.0,
            "lower_bound" => OP.HybridDualDynamicProgrammingAlgo(
                QP_SOLVER,
                CDDLib.Library(),
                1e-5,
                1e-4,
                1,
            ),
        ),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), _glb_problem())
    MOI.optimize!(optimizer)
    return optimizer
end

# UniformEllipsoidAbstraction on the PWA system (SDP transitions via Clarabel, two phases).
function run_uniform_ellipsoid_pwa()
    Random.seed!(0)
    lib = CDDLib.Library()
    dt = 0.01
    concrete_problem =
        PwaSystem.problem(; lib = lib, dt = dt, Usz = 70, Wsz = 3, simple = true)
    concrete_system = concrete_problem.system

    alternating_simulation_problem =
        PR.AlternatingSimulationProblem(concrete_system, concrete_system.ext[:X])

    n_step = 3
    origin = SVector(0.0, 0.0)
    h = SVector(1.0 / n_step, 1.0 / n_step)
    nx = size(concrete_system.resetmaps[1].A, 1)
    P = (1 / nx) * diagm((h ./ 2) .^ (-2))
    Pm = P
    R = h ./ 2
    state_grid = MP.GridEllipsoidalRectangular(origin, h, P)

    optimizer = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("alternating_simulation_problem"),
        alternating_simulation_problem,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), MP.INNER)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("P"), P)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("Pm"), Pm)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("R"), R)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), SDP_SOLVER)
    naug = nx + 2 + 1  # nx + nu + 1, nu = 2 for this instance
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("Q_aug"),
        Matrix{Float64}(I, naug, naug) * (dt^2),
    )

    MOI.optimize!(optimizer)  # phase 1: build abstraction
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.optimize!(optimizer)  # phase 2: synthesize controller
    return optimizer
end

# LazyEllipsoidsAbstraction on the nonlinear system (RRT + SDP/Lyapunov via Clarabel).
function run_lazy_ellipsoids_nonlinear()
    Random.seed!(0)
    concrete_problem = NonLinear.problem()
    obstacles = NonLinear.default_obstacles()

    optimizer = MOI.instantiate(AB.LazyEllipsoidsAbstraction.Optimizer)
    AB.LazyEllipsoidsAbstraction.set_optimizer!(
        optimizer,
        concrete_problem,
        SDP_SOLVER,
        100,      # maxδx
        10 * 2,   # maxδu
        0.01,     # λ
        1,        # k1
        1,        # k2
        false,    # RRTstar
        false,    # continues
        100;      # maxIter
        obstacles = obstacles,
    )
    MOI.optimize!(optimizer)
    return optimizer
end

# ---------------------------------------------------------------------------
# The suite: (solver family, problem) pairs, ordered fast-first.
# ---------------------------------------------------------------------------
const CASES = [
    (
        method = "UniformGridAbstraction",
        problem = "DCDC (safety)",
        run = run_uniform_grid_dcdc,
    ),
    (
        method = "UniformGridAbstraction",
        problem = "PathPlanning (reach-avoid)",
        run = run_uniform_grid_path_planning,
    ),
    (
        method = "BemporadMorari",
        problem = "GolLazarBelta (N=9)",
        run = run_bemporad_morari_glb,
    ),
    (
        method = "BranchAndBound",
        problem = "GolLazarBelta (N=9)",
        run = run_branch_and_bound_glb,
    ),
    (
        method = "UniformEllipsoidAbstraction",
        problem = "PwaSystem",
        run = run_uniform_ellipsoid_pwa,
    ),
    (
        method = "LazyEllipsoidsAbstraction",
        problem = "NonLinear",
        run = run_lazy_ellipsoids_nonlinear,
    ),
]
