using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
using LinearAlgebra
using JLD2
import MathOptInterface as MOI

include("../robot_vcontrol.jl")
import .RobotVelocity as RV

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems


# ------------------------------------------------------------
# 1) Define a simple 4D continuous-time system: x' = u
# ------------------------------------------------------------

robot_geometry = RV.Geometry.geometry(0.0, 0.202, 0.172)

# ------------------------------------------------------------
# 2) Abstraction construction (AlternatingSimulationProblem)
# ------------------------------------------------------------

# State space
x_bar = 1.5
full_X_ = UT.HyperRectangle(
    -SVector(x_bar, x_bar, x_bar, x_bar),
    SVector(x_bar, x_bar, x_bar, x_bar)
)

# State grid resolution (hx defines TOTAL CELL WIDTH)
dx = 0.1
x0 = SVector(0.0, 0.0, 0.0, 0.0)
hx = SVector(dx, dx, dx, dx)
state_grid = MP.GridFree(x0, hx)

obstacle = UT.HyperRectangle(SVector(-0.03, 0.0), SVector(0.03,0.03))

# State constraints
_X_ = RV.RobotModel.remove_infeasible_cells(
    full_X_,
    state_grid,
    obstacle,
    robot_geometry, true
)
# _X_ = full_X_

# Input space
u_bar = 10*dx
_U_ = UT.HyperRectangle(
    -SVector(u_bar, u_bar, u_bar, u_bar),
    SVector(u_bar, u_bar, u_bar, u_bar)
)

# Input grid resolution
du = 10*dx
u0 = SVector(0.0, 0.0, 0.0, 0.0)
hu = SVector(du, du, du, du)
input_grid = MP.GridFree(u0, hu)

# Initial set
_I_ = UT.HyperRectangle(
    SVector(0.2, 0.0, -0.2, 0.0) - hx/2,
    SVector(0.2, 0.0, -0.2, 0.0) + hx/2
)

_T_ = RV.RobotModel.compute_target_set(state_grid, SVector(0.25, 0.0), robot_geometry, true)

concrete_system = RV.RobotModel.system(; _X_ = _X_, _U_ = _U_)

# concrete_problem = DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

# ------------------------------------------------------------
# 3) Define co-safe LTL problem with sets labeling
# ------------------------------------------------------------

concrete_problem =
    DI.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        0.0
    )