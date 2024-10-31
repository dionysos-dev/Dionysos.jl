using JuMP
using StaticArrays
using CDDLib
using LinearAlgebra
using Clarabel
using Mosek, MosekTools
using OSQP
using Ipopt
using HiGHS
using Pavito

using Dionysos
const OP = Dionysos.Optim
const AB = OP.Abstraction
const PB = Dionysos.Problem
const DO = Dionysos.Domain
const UT = Dionysos.Utils
const SY = Dionysos.Symbolic

bench = Dict{Tuple{String, String}, MOI.AbstractOptimizer}()

#######################
# Import all problems #
#######################

PROBLEMS_PATH = joinpath(@__DIR__, "..", "problems")
PROBLEMS_AVOID = ["robot_problem.jl"]
problems_modules = Dict{String, Module}()
for file in readdir(PROBLEMS_PATH)
    if file ∈ PROBLEMS_AVOID || !endswith(file, ".jl")
        continue
    end
    path = joinpath(PROBLEMS_PATH, file)
    mod = include(path)
    problems_modules[file] = mod
end

###################
# Bemporad Morari #
###################

osqp = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true,
)

ipopt = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true);
QP_SOLVERS = [("OSQP", osqp), ("Ipopt", ipopt)]

highs = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true);
MIP_SOLVERS = [("HiGHS", highs)]

MIQP_SOLVERS = []
for (mip_name, mip_solver) in MIP_SOLVERS
    for (qp_name, qp_solver) in QP_SOLVERS
        pavito = optimizer_with_attributes(
            Pavito.Optimizer,
            "mip_solver" => mip_solver,
            "cont_solver" => qp_solver,
            MOI.Silent() => true,
        )
        push!(MIQP_SOLVERS, ("Pavito($mip_name,$qp_name)", pavito))
    end
end

for (qp_name, qp_solver) in QP_SOLVERS
    for (miqp_name, miqp_solver) in MIQP_SOLVERS
        solver_bm = optimizer_with_attributes(
            OP.BemporadMorari.Optimizer{Float64},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0,
            "problem" => problems_modules["gol_lazar_belta.jl"].problem(),
        )
        bench["BemporadMorari($qp_name, $miqp_name)", "gol_lazar_belta.jl"] =
            MOI.instantiate(solver_bm)
    end
end

####################
# UniformGridAbstraction #
####################

x0 = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)
DCDC = problems_modules["dc_dc.jl"]
solver_uniform_grid_dc_dc = optimizer_with_attributes(
    AB.UniformGridAbstraction.Optimizer,
    "concrete_problem" => DCDC.problem(),
    "state_grid" => state_grid,
    "input_grid" => input_grid,
    "jacobian_bound" => DCDC.jacobian_bound(),
    "time_step" => 0.5,
)
bench["UniformGridAbstraction", "dc_dc.jl"] = MOI.instantiate(solver_uniform_grid_dc_dc)

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
state_grid = DO.GridFree(x0, h);
u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, h);
PathPlanning = problems_modules["path_planning.jl"]
solver_uniform_grid_path_planning = optimizer_with_attributes(
    AB.UniformGridAbstraction.Optimizer,
    "concrete_problem" => PathPlanning.problem(; simple = true),
    "state_grid" => state_grid,
    "input_grid" => input_grid,
    "jacobian_bound" => PathPlanning.jacobian_bound(),
    "time_step" => 0.3,
)
bench["UniformGridAbstraction", "path_planning.jl"] =
    MOI.instantiate(solver_uniform_grid_path_planning)

###################
# LazyAbstraction #
###################

# specific functions
function post_image(abstract_system, concrete_system, xpos, u)
    Xdom = abstract_system.Xdom
    x = DO.get_coord_by_pos(Xdom.grid, xpos)
    Fx = concrete_system.f_eval(x, u)
    r = Xdom.grid.h / 2.0 + concrete_system.measnoise
    Fr = r

    rectI = DO.get_pos_lims_outer(Xdom.grid, UT.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(DO._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        ypos = DO.set_in_period_pos(Xdom, ypos)
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = SY.get_state_by_xpos(abstract_system, ypos)[1]
        push!(over_approx, target)
    end
    return allin ? over_approx : []
end

function pre_image(abstract_system, concrete_system, xpos, u)
    grid = abstract_system.Xdom.grid
    x = DO.get_coord_by_pos(grid, xpos)
    potential = Int[]
    x_prev = concrete_system.f_backward(x, u)
    xpos_cell = DO.get_pos_by_coord(grid, x_prev)
    n = 2
    for i in (-n):n
        for j in (-n):n
            x_n = (xpos_cell[1] + i, xpos_cell[2] + j)
            x_n = DO.set_in_period_pos(abstract_system.Xdom, x_n)
            if x_n in abstract_system.Xdom
                cell = SY.get_state_by_xpos(abstract_system, x_n)[1]
                if !(cell in potential)
                    push!(potential, cell)
                end
            end
        end
    end
    return potential
end

function compute_reachable_set(rect::UT.HyperRectangle, concrete_system, Udom)
    r = (rect.ub - rect.lb) / 2.0 + concrete_system.measnoise
    Fr = r
    x = UT.get_center(rect)
    n = UT.get_dims(rect)
    lb = fill(Inf, n)
    ub = fill(-Inf, n)
    for upos in DO.enum_pos(Udom)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        Fx = concrete_system.f_eval(x, u)
        lb = min.(lb, Fx .- Fr)
        ub = max.(ub, Fx .+ Fr)
    end
    lb = SVector{n}(lb)
    ub = SVector{n}(ub)
    return UT.HyperRectangle(lb, ub)
end
minimum_transition_cost(symmodel, contsys, source, target) = 1.0

problem_simple = problems_modules["simple_problem.jl"].problem()
concrete_system = problem_simple.system

hx = [0.5, 0.5]
u0 = SVector(0.0, 0.0)
hu = SVector(0.5, 0.5)
Ugrid = DO.GridFree(u0, hu)
hx_heuristic = [1.0, 1.0] * 1.5
maxIter = 100

solver_lazy_simple = MOI.instantiate(AB.LazyAbstraction.Optimizer)
AB.LazyAbstraction.set_optimizer!(
    solver_lazy_simple,
    problem_simple,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    minimum_transition_cost,
    hx_heuristic,
    hx,
    Ugrid,
)

bench["LazyAbstraction", "simple_problem.jl"] = solver_lazy_simple

###################
# LazyAbstraction #
###################

problem_simple = problems_modules["simple_problem.jl"].problem(;
    rectX = UT.HyperRectangle(SVector(0.0, 0.0), SVector(60.0, 60.0)),
    obstacles = [UT.HyperRectangle(SVector(22.0, 21.0), SVector(25.0, 32.0))],
    periodic = Int[],
    periods = [30.0, 30.0],
    T0 = [0.0, 0.0],
    rectU = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0)),
    Uobstacles = [UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))],
    _I_ = UT.HyperRectangle(SVector(6.5, 6.5), SVector(7.5, 7.5)),
    _T_ = UT.HyperRectangle(SVector(44.0, 43.0), SVector(49.0, 48.0)),
    state_cost = UT.ZeroFunction(),
    transition_cost = UT.ConstantControlFunction(1.0),
    tstep = 0.8,
    measnoise = SVector(0.0, 0.0),
)
concrete_system = problem_simple.system

hx_local = [0.5, 0.5]
hx_heuristic = [1.0, 1.0]
u0 = SVector(0.0, 0.0)
hu = SVector(0.5, 0.5)
Ugrid = DO.GridFree(u0, hu)

local_optimizer = MOI.instantiate(AB.LazyAbstraction.Optimizer)
AB.LazyAbstraction.set_optimizer_parameters!(
    local_optimizer,
    100,
    pre_image,
    post_image,
    compute_reachable_set,
    minimum_transition_cost,
    hx_local,
    hx_heuristic;
    γ = 10.0,
)

hx_global = [10.0, 10.0]
u0 = SVector(0.0, 0.0)
hu = SVector(0.5, 0.5)
Ugrid = DO.GridFree(u0, hu)
max_iter = 6
max_time = 1000

solver_hierarchical_simple = MOI.instantiate(AB.HierarchicalAbstraction.Optimizer)
AB.HierarchicalAbstraction.set_optimizer!(
    solver_hierarchical_simple,
    problem_simple,
    hx_global,
    Ugrid,
    compute_reachable_set,
    minimum_transition_cost,
    local_optimizer,
    max_iter,
    max_time;
    option = true,
)

bench["HierarchicalAbstraction", "simple_problem.jl"] = solver_hierarchical_simple

#########################
# EllipsoidsAbstraction #
#########################

lib = CDDLib.Library() # polyhedron lib
Usz = 70 # upper limit on |u|, `Usz = 50` in [1]
Wsz = 3 # `Wsz = 5` in [1]
dt = 0.01; # discretization step

problem_pwa = problems_modules["pwa_sys.jl"].problem(;
    lib = lib,
    dt = dt,
    Usz = Usz,
    Wsz = Wsz,
    simple = false,
)
concrete_system = problem_pwa.system

n_step = 3
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0 / n_step, 1.0 / n_step)
nx = size(concrete_system.resetmaps[1].A, 1)
P = (1 / nx) * diagm((X_step ./ 2) .^ (-2))
state_grid = DO.GridEllipsoidalRectangular(X_origin, X_step, P, concrete_system.ext[:X]);
opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

solver_ellispoids_pwa = optimizer_with_attributes(
    AB.EllipsoidsAbstraction.Optimizer,
    "concrete_problem" => problem_pwa,
    "state_grid" => state_grid,
    "sdp_solver" => opt_sdp,
)

bench["EllipsoidsAbstraction", "pwa_sys.jl"] = MOI.instantiate(solver_ellispoids_pwa)

#############################
# LazyEllipsoidsAbstraction #
#############################

problem_nonlinear = problems_modules["non_linear.jl"].problem()
concrete_system = problem_nonlinear.system

# Optimizer's parameters
const FALLBACK_URL = "mosek://solve.mosek.com:30080"
sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
MOI.set(sdp_opt, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)

maxδx = 100
maxδu = 10 * 2
λ = 0.01
k1 = 1
k2 = 1
RRTstar = false
continues = false
maxIter = 100

solver_lazy_ellipsoids_nonlinear = MOI.instantiate(AB.LazyEllipsoidsAbstraction.Optimizer)
AB.LazyEllipsoidsAbstraction.set_optimizer!(
    solver_lazy_ellipsoids_nonlinear,
    problem_nonlinear,
    sdp_opt,
    maxδx,
    maxδu,
    λ,
    k1,
    k2,
    RRTstar,
    continues,
    maxIter,
)

bench["LazyEllipsoidsAbstraction", "non_linear.jl"] = solver_lazy_ellipsoids_nonlinear
