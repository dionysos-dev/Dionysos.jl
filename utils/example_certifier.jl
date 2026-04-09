using StaticArrays
using JuMP
using Plots
import MathOptInterface as MOI

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction
const SC = AB.SymbolicCertifier

include("../problems/toy_problem.jl")

# ------------------------------------------------------------
# 1) Build the optimizer
# ------------------------------------------------------------
concrete_problem = ToyProblem.optimal_control_problem()

# State discretization
hx = SVector(0.1, 0.1)
# Input discretization
u0 = SVector(0.0, 0.0)
hu = SVector(0.2, 0.2)
# Time discretization
tstep = 0.2

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), ToyProblem.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> ST.NewIndexedAutomatonList(n, m),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

# ------------------------------------------------------------
# 2) Build a *candidate* trajectory (simple heuristic)
# ------------------------------------------------------------
candidate_x_traj = Dionysos.System.Trajectory{SVector{2, Float64}}(
    SVector{2, Float64}[
        [-1.5, -1.5],
        [-1.5, -1.5],
        [-0.31032202176824586, -0.31032202176824586],
        [-0.2452650817125079, -0.2452650817125079],
        [-0.13999059290084434, -0.13999059290084434],
    ],
)

# ------------------------------------------------------------
# 3) Call your local tube certifier
# ------------------------------------------------------------
cert = SC.UniformGridLocalTubeCertifier()

SC.set_optimizer!(cert, optimizer)
SC.set_trajectory!(cert, candidate_x_traj)
cert.radius = SVector(0.25, 0.25)
cert.margin = 0.0
cert.incl_mode = MP.OUTER
cert.enforce_safe_max_step = true

SC.certify!(cert)

println("\n=== Local Certification Result ===")
println("success:    ", SC.get_success(cert))
println("time (sec): ", SC.get_solve_time(cert))
controller = SC.get_controller(cert)

# ------------------------------------------------------------
# 4) Closed loop trajectory 
# ------------------------------------------------------------

x0 = SVector(UT.get_center(concrete_problem.initial_set)...)
reached(x) = (x ∈ concrete_problem.target_set)
nstep = 300
x_traj, u_traj = ST.get_closed_loop_trajectory(
    MOI.get(cert.optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    controller,
    x0,
    nstep;
    stopping = reached,
);

# ------------------------------------------------------------
# 5) Plot 
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal, title = "ToyExample: candidate traj + sets")
plot!(concrete_problem.system.X; color = :grey, opacity = 0.15, label = "X")
plot!(concrete_problem.initial_set; color = :green, opacity = 0.25, label = "Initial set")
plot!(concrete_problem.target_set; color = :red, opacity = 0.35, label = "Target set")
plot!(
    cert.optimizer.abstraction_solver.abstract_system;
    color = :blue,
    opacity = 0.35,
    label = "Tube",
    efficient = false,
)
plot!(candidate_x_traj; ms = 2.0, arrows = false, label = "Candidate")
plot!(x_traj; color = :blue, ms = 2.0, arrows = false, label = "Closed loop Trajecory")
display(fig)
