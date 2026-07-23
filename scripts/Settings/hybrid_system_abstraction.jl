using StaticArrays: SVector, SMatrix
import MathOptInterface as MOI

using Plots

import MathematicalSystems as MS
using HybridSystems

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

# ==============================
# Define the HybridSystem
# ==============================

# State and input sets (1D)
X = UT.box(SVector(-1.0), SVector(1.0))
U = UT.box(SVector(-1.5), SVector(1.5))

# Mode dynamics (return SVector for StaticArrays friendliness)
mode1_f(x, u) = SVector(0.5 * x[1] + u[1])
mode2_f(x, u) = SVector(0.8 * x[1] + u[1])

mode1_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
mode2_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

# Time system (1D time in [0,3], derivative 1)
Tdom = UT.box(SVector(0.0), SVector(3.0))
time_sys = MS.ConstrainedLinearContinuousSystem(SMatrix{1, 1}(1.0), Tdom)

# Guard and reset for switching (augmented state is [x; t] => 2D)
guard_1 = UT.box(SVector(0.2, 0.0), SVector(0.7, 0.9))

struct FixedPointResetMap <: MS.AbstractMap
    domain::UT.Box
    target::Vector{Float64}
end

MS.apply(reset::FixedPointResetMap, state::AbstractVector) = [reset.target[1], state[2]]
MS.stateset(reset::FixedPointResetMap) = reset.domain

reset_map = FixedPointResetMap(guard_1, SVector(0.0, 0.0))  # reset to (x,t)=(0, t)

# Automaton + hybrid system
automaton = HybridSystems.GraphAutomaton(2)
HybridSystems.add_transition!(automaton, 1, 2, 1)

modes_systems = [
    ST.VectorContinuousSystem([mode1_system, time_sys]),
    ST.VectorContinuousSystem([mode2_system, time_sys]),
]
reset_maps = [reset_map]
switchings = [HybridSystems.AutonomousSwitching()]

concrete_system =
    HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

# ==============================
# Define the Control Problem
# ==============================

# Augmented initial state: (x, t, mode)
initial_state = (SVector(0.0), 0.0, 1)

# Target: in mode 2, x ∈ [-1,1], t ∈ [1,2]
target_set = PR.hybrid_reach_spec(
    [UT.box(SVector(-1.0), SVector(1.0))],
    [UT.box(SVector(1.0), SVector(2.0))],
    [2],
)

trans_cost = (aug_state, u) -> 1.0

concrete_problem = PR.OptimalControlProblem(
    concrete_system,
    initial_state,
    target_set,
    nothing,
    trans_cost,
)

# ==============================
# Solver parameters
# ==============================

# Per-mode abstraction subsolvers: pass *constructors* (zero-arg callables)
optimizer_list = [
    () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
]

# Grids (GridFree signature depends on your Mapping API; these are SVectors)
state_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
input_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
state_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))
input_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))

optimizer_kwargs_dict = [
    Dict{String, Any}(
        "state_grid" => state_grid_1,
        "input_grid" => input_grid_1,
        "time_step" => 0.1,
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => (u) -> SMatrix{1, 1, Float64}(0.5),
    ),
    Dict{String, Any}(
        "state_grid" => state_grid_2,
        "input_grid" => input_grid_2,
        "time_step" => 0.1,
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => (u) -> SMatrix{1, 1, Float64}(0.8),
    ),
]

# Main solver
optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
    optimizer_kwargs_dict,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

# Solve
MOI.optimize!(optimizer)

# Retrieve results
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# ==============================
# Closed-loop Trajectory
# ==============================

reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

max_steps = 20
tsteps = [0.1, 0.1] # per mode
aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    tsteps,
    initial_state,
    max_steps;
    stopping = reached,
)

# aug_state = (x, t, mode)
tx_mode1 = ST.Trajectory([SVector(t, x[1]) for (x, t, k) in aug_x_traj if k == 1])
tx_mode2 = ST.Trajectory([SVector(t, x[1]) for (x, t, k) in aug_x_traj if k == 2])

fig = plot(; aspect_ratio = :equal)
plot!(fig, tx_mode1; dims = [1, 2], label = "mode 1", color = :blue)
plot!(fig, tx_mode2; dims = [1, 2], label = "mode 2", color = :red)
display(fig)
