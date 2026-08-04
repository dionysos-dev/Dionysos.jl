# Sizing a center-simulation uniform-grid abstraction from a δ-GAS Lyapunov function.
#
# In the classical uniform-grid abstraction with `CENTER_SIMULATION`, only the center of
# each cell is propagated. That is unsound on its own — but for an incrementally stable
# (δ-GAS) system, a common quadratic Lyapunov function bounds the intra-cell drift, which
# both *sets the grid step* and gives the precision ε of the resulting
# ε-approximately-bisimilar symbolic model (Girard–Pola–Tabuada, IEEE TAC 2010).
#
# The DC-DC boost converter is the canonical example: two switched affine modes, the input
# selecting the mode. Here we compute the δ-GAS certificate, derive the grid from a target
# precision, then build the center-simulation abstraction and synthesize a safety
# controller.

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using JuMP
import MathOptInterface as MOI
import StaticArrays: SVector
import Clarabel
import LazySets
import HybridSystems
const PCLF = UT.PathCompleteFramework

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC", "dcdc_converter.jl"),
)

# --- 1. Common quadratic (δ-GAS) Lyapunov function of the sampled modes ---------------
# The error e = x − y obeys the same switched dynamics, so a common quadratic Lyapunov
# function for the *sampled* modes exp(A_σ·τ) is exactly the δ-GAS certificate. We reuse
# the path-complete quadratic SDP with a single-node (common-Lyapunov) graph.
params = DCDC.Params()
τ = 0.5                                    # sampling time
sampled_modes = [exp(Matrix(DCDC.A1(params)) * τ), exp(Matrix(DCDC.A2(params)) * τ)]
switched_system = HybridSystems.discreteswitchedsystem(sampled_modes)

sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
pclf = PCLF.compute_quadratic_pieces_pclf(
    switched_system,
    PCLF.generate_DeBruijn_edges(length(sampled_modes), 0),   # single-node (common) graph
    sdp_optimizer;
    MLF = true,
)
P = pclf.pieces[1].P            # the single common Lyapunov matrix
γ = pclf.JSRapprox              # per-step contraction of V(x, y) = ‖x - y‖_P
println(
    "Per-step contraction γ = ",
    round(γ; sigdigits = 4),
    "  (equivalent rate κ = -log(γ)/τ = ",
    round(-log(γ) / τ; sigdigits = 4),
    ")",
)
println("Lyapunov matrix P = ", round.(P; sigdigits = 4))

# --- 2. Grid step ↔ precision -------------------------------------------------------
concrete_problem = DCDC.problem()
X = concrete_problem.system.X
widths = LazySets.high(X) .- LazySets.low(X)

# GPT sizing: what grid does a target precision demand, and how many cells is that?
println("\nprecision ε   →   grid step η   →   ≈ #cells over X")
for ε in (0.05, 0.1, 0.2)
    η = UT.grid_step_from_lyapunov(P, ε, γ)
    ncells = prod(ceil.(Int, widths ./ η))
    println("  ε = $ε   →   η = $(round(η; sigdigits = 3))   →   ≈ $ncells cells")
end

# The classical hand-picked grid hx = 2/4000 corresponds to a definite precision:
η_classic = 2.0 / 4.0e3
println(
    "\nClassical grid η = $η_classic certifies precision ε = ",
    round(UT.precision_from_grid_step(P, η_classic, γ); sigdigits = 4),
)

# For a quick, tractable demo we use a coarser grid (fewer cells) and report the
# precision it certifies. Shrink η toward η_classic for a tighter (but heavier) model.
η = 5.0e-3
ε_demo = UT.precision_from_grid_step(P, η, γ)
println(
    "\nDemo grid η = $η  ⇒  the center-simulation abstraction is ",
    "ε = $(round(ε_demo; sigdigits = 4))-approximately bisimilar.",
)

# --- 3. Build the center-simulation abstraction with that η and synthesize safety -----
state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(η, η))
input_grid = MP.GridFree(SVector(1), SVector(1))          # inputs {1, 2} select the mode

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), τ)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

println(
    "\nAbstraction: $(SY.get_n_state(abstract_system)) states, ",
    "$(SY.get_n_input(abstract_system)) inputs.",
)
println(
    "Synthesized safety controller is ε = $(round(ε_demo; sigdigits = 4))-approximately ",
    "correct on the concrete system (Girard–Pola–Tabuada).",
)

# --- 4. Closed-loop sanity check ------------------------------------------------------
x0 = SVector(1.2, 5.6)
traj = ST.get_closed_loop_trajectory(discrete_time_system, concrete_controller, x0, 100)
println("Closed-loop trajectory length: ", length(collect(ST.states(traj))), " steps.")

# --- 5. Dashboard visualization of the closed-loop safety run -------------------------
system_plot! = DCDC.system_plot!()
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),
    udims = (1,),
    Δt = τ,
    fps = 5,
    # filename = "dcdc_incremental_dashboard.mp4",
    xlabel_state = "iL",
    ylabel_state = "vC",
    xlabel_input = "time [s]",
    ylabel_input = "mode",
    xlims_state = (1.15, 1.55),
    ylims_state = (5.45, 5.85),
    ylims_input = (0.5, 2.5),
)
