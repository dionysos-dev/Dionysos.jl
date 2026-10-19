# # Underactuated four-link biped in normal form: a symmetric forward step
#
# One footstep of a four-link planar biped, with the three motor joints commanded in velocity and the
# ankle `q₁` unactuated. **The default experiment is a symmetric forward step**: the swing foot starts at
# `x ≈ -0.2 m` behind the stance foot and has to land around the foothold at `x = +0.2 m` in front of it,
# with the ankle falling forwards all along (the ankle momentum stays non-positive and ends at any
# negative value), without obstacle for now (`OBSTACLE = :bump` brings the bump back). The model and the
# abstraction are:
#
# - **model**: the multibody system is built *symbolically* by `RP.symbolic_model()`
#   (`problems/BipedRobot/4D_model_underactuated/robot_problem.jl`, on top of `Multibody2D.jl`; no
#   RigidBodyDynamics in the dynamics), brought to normal form `x = (q, v_a, p_u)`, and reduced to the
#   velocity-controlled design model in which the unactuated ankle is described by its
#   **generalized momentum** `p₁ = ∂L/∂q̇₁` (`RP.reduced_dynamics`):
#
#       x = [p₁, q₁, q₂, q₃, q₄],   u = [v₂, v₃, v₄],
#
#   i.e. `q̇₁ = (p₁ - D₁ₐ v_a) / D₁₁`, `ṗ₁ = ∂L/∂q₁`, `(q̇₂, q̇₃, q̇₄) = v_a`, with the commands
#   held constant over a sampling period;
# - **abstraction**: a certified over-approximation with the **local growth bound** of
#   `LocalGrowthAbstraction.jl` (a Jacobian bound per cell and input, with interval arithmetic and
#   a verified flow enclosure). The abstraction is nondeterministic, so the controller is synthesized
#   against the worst case.
#
# The constraints: the swing-foot height is bounded below by the ground (cells where the foot is
# surely below `z = 0` are removed, the constraint named `"ground"`), the ankle momentum is bounded
# above by `0` (a forward step: the state box stops there), and, with `OBSTACLE = :bump`, a triangular
# obstacle at the stance foot is carved out, which forces the foot to lift (`"lift constraint"`).
#
# Other experiments are drivers of this file: a variable defined before the `include` replaces the
# default of the same name.
#
# Conventions: `Multibody2D` and `robot_problem.jl` both measure angles counter-clockwise, so the
# symbolic model is used as is. The forward kinematics of the carving and of the target set are
# compiled from the same symbolic model (`RP.Geometry`), so the dynamics and the constraints cannot
# disagree. The mechanism is the one of `four_link_biped.urdf` (thin rods `m L² / 12`, point-mass
# motors).
#
# **Cost.** The abstraction has one image to compute per `(cell, input)` pair (27 inputs by default), and the
# number of cells is the *product* of the cell counts along the five axes: a state box that looks only
# slightly wider can have ten times more cells (the size is printed before the abstraction starts).
#
# Run it with (an environment with Dionysos, Symbolics, IntervalArithmetic, Plots and the other
# packages of `problems/BipedRobot`; RigidBodyDynamics is not needed):
#
# ```
# julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_normal_form.jl
# ```

using StaticArrays
using MathematicalSystems
using Dionysos
import MathOptInterface as MOI
import LazySets

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

const PROBLEM_DIR = joinpath(
    dirname(dirname(pathof(Dionysos))),
    "problems",
    "BipedRobot",
    "4D_model_underactuated",
)

# Mechanism, geometry, obstacle carving, target set and plots (state `[ankle, q₁, q₂, q₃, q₄]`).
include(joinpath(PROBLEM_DIR, "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# Local growth bound.
include(joinpath(PROBLEM_DIR, "LocalGrowthAbstraction.jl"))
import .LocalGrowthAbstraction as LGA

# ------------------------------------------------------------
# 1) Symbolic model → numerical dynamics
# ------------------------------------------------------------

values = RP.default_parameters()   # lengths, masses and gravity of four_link_biped.urdf

@time model = RP.symbolic_model()
println("Symbolic normal form: ", model.nf)
reduced = RP.reduced_dynamics(model)       # x = (p₁, q₁…q₄), u = (v₂, v₃, v₄), parameters still symbolic
f_sym = RP.substitute_parameters(model, reduced.f, values)   # now only x and u are free
x_sym, u_sym = reduced.x, reduced.u

# The forward kinematics of the carving and of the target set, compiled from the same model.
geometry = RP.Geometry(model, values)

# ------------------------------------------------------------
# 2) Scenario: a symmetric forward step
# ------------------------------------------------------------
# A foothold 20 cm in front of the stance foot, and (optionally) an 8 cm × 3 cm bump at the stance
# foot that the swing foot must clear.

_triangle(cx, hw, h) =
    LazySets.VPolygon([SVector(cx - hw, 0.0), SVector(cx, h), SVector(cx + hw, 0.0)])

# `SCENARIO_NAME` only changes the names of the output files. `OBSTACLE` is `nothing` (none, the
# default), `:bump` (the triangle above) or a vector of `LazySets` polygons.
scenario_name = @isdefined(SCENARIO_NAME) ? SCENARIO_NAME : :symmetric_step
obstacle = @isdefined(OBSTACLE) ? OBSTACLE : nothing
obstacle === :bump && (obstacle = [_triangle(0.0, 0.04, 0.03)])
foothold = SVector(0.2, 0.0)

# The swing foot must land within `FOOTHOLD_TOLERANCE` metres of the foothold along `x` (and up to
# that height above the ground). With `0.0` it has to land *exactly* on it, which is a thin sliver of
# the angle space: the over-approximated ankle spreads over two or three cells every step, and a
# controller certified against the worst case cannot land inside a sliver, whatever the start. A few
# centimetres of tolerance make the problem solvable: 3 cm for a step of 19 cm, **8 cm** for the 40 cm
# step of this experiment on this coarse grid (the foot lands anywhere in `x ∈ [0.12, 0.28]`, up to 8 cm above the ground).
foothold_tolerance = @isdefined(FOOTHOLD_TOLERANCE) ? FOOTHOLD_TOLERANCE : 0.05

# Set `DISC` before including this file to try another resolution. The motor joints move by
# exactly one cell (`dx`) per step at speed `du = dx / tstep`; `dq̇1` is here the step of the
# ankle *momentum* `p₁` (`D₁₁ ≈ 0.14 kg m²`, so 0.04 corresponds to ≈ 0.3 rad/s). The default is a
# **short period** (`tstep = 1/60 s`, hence motors at `du = 6 rad/s`) on a **coarse** grid
# (`dq1 = 0.1`, `dq̇1 = 0.04`) to keep the run short while the box is large: every step the
# over-approximated ankle spreads over several cells whatever the motion, so the controller needs
# frequent corrections. Refine the grid (`dq1 = 0.05`, `dq̇1 = 0.02`) once the start is chosen: the
# cells are 4 times more numerous and the required tolerance smaller.
#
# `input_steps = k` is the number of speeds each motor can take on each side: the commands are the integer
# multiples of `du` in `{-k du, …, -du, 0, du, …, k du}` (`2k + 1` values per joint, `(2k + 1)³` inputs: 27
# for `k = 1`, 125 for `k = 2`, 343 for `k = 3`), so a joint moves by up to `k` cells per step. The cost of the
# abstraction is proportional to the number of inputs, and the ankle momentum changes more per step
# with the larger commands.
disc =
    @isdefined(DISC) ? DISC :
    RP.default_discretization(;
        dx = 0.05,
        dq1 = 0.05,
        dq̇1 = 0.04,
        tstep = 1 / 60,
        input_steps = 3,
    )

# State box `[p₁, q₁, q₂, q₃, q₄]`. **Its size is the cost of the run**: the number of cells is the
# product of the cell counts along the five axes (this box: 123 420 cells, 3.3 M cell/input pairs:
# about 2 minutes on 10 threads), and a box that looks only slightly wider can have ten times more (the
# size is printed before the abstraction). The step is **forward**: the ankle momentum is bounded
# above by `-0.08`, which **forces the ankle to keep a forward momentum along the whole step** (a
# transition that would leave the box is forbidden), and below by `-0.6` (`p₁ ≈ D₁₁ q̇₁ + D₁ₐ v_a`,
# `D₁₁ ≈ 0.14`, so `q̇₁ ≳ -4 rad/s`). The step starts from a large negative momentum: with the box
# stopping at `-0.2` the momentum goes back towards 0 before the end of the step (the robot starts
# falling back), and cells too close to the bounds have no controller because the over-approximation
# can push them out. Enlarging the angle bounds (e.g. `q₁ ∈ [-0.5, 0.3]`, `q₂ ∈ [-0.8, 0.8]`,
# `q₃ ∈ [-1, 1.2]`, `q₄ ∈ [-1.2, 0.4]`) gives more starts but 4 times more cells. Angles are
# counter-clockwise: the swing hip swings forwards for `q₃ > 0`, the swing knee bends backwards for
# `q₄ < 0`, and the robot falls forwards, towards negative `q₁`. Set `X_BOUNDS = (x_lb, x_ub)` to
# change it.
x_lb, x_ub =
    @isdefined(X_BOUNDS) ? X_BOUNDS :
    (SVector(-0.6, -0.4, -1.0, -1.0, -1.0), SVector(-0.08, 0.2, 0.1, 1.0, 0.1))

# The footstep starts from the configuration `q0 = [q₁, q₂, q₃, q₄]` and ends on the foothold. The
# controller must work from every cell of the ankle momentum `p₁` that `INITIAL_VELOCITY` touches and
# may end with any value in `TARGET_VELOCITY` (`nothing`, the default: the whole interval of
# definition, here the negative momenta of the box). Both are `(vlo, vhi)` ranges on `p₁`. The default
# starts from **a single value**, `p₁ = -0.4` (one cell of width 0.04: the ankle is already falling
# forwards fast, at about -2.5 rad/s), with the swing foot at `(-0.2, 0.028)` m (stance femur tilted
# 0.5 rad backwards, swing femur vertical, swing knee bent backwards by 0.9 rad): the step is symmetric,
# from `x ≈ -0.2` to the foothold at `+0.2`. Set `Q0`, `INITIAL_VELOCITY` or `TARGET_VELOCITY` before
# including this file to change them, e.g. `INITIAL_VELOCITY = (-0.44, -0.36)`.
#
# With an obstacle, the swing foot has to start **outside** its carved neighbourhood: the carving is
# sound, so it removes every cell whose swing foot is within the Lipschitz margin (printed below) of
# the obstacle.
q0 = @isdefined(Q0) ? Q0 : SVector(0.0, 0.5, -0.4, -0.9)
initial_velocity = @isdefined(INITIAL_VELOCITY) ? INITIAL_VELOCITY : (-0.4, -0.4)
target_velocity = @isdefined(TARGET_VELOCITY) ? TARGET_VELOCITY : nothing

# All combinations of the three velocity commands by default (27 for `input_steps = 1`): keeping the ankle
# balanced while the swing leg moves needs coordinated joints. Set `ONE_JOINT_PER_STEP = true` before
# including this file to allow one moving motor joint per step (7 inputs for `input_steps = 1`, a smaller automaton).
one_joint_per_step = @isdefined(ONE_JOINT_PER_STEP) ? ONE_JOINT_PER_STEP : false
extra_attributes =
    one_joint_per_step ?
    (; state_input_filter = (x, u) -> count(v -> abs(v) > 1e-9, u) <= 1) : (;)

println("Carving the obstacle and the ground out of the state domain…")
X_box = LazySets.Hyperrectangle(; low = x_lb, high = x_ub)
U_box = LazySets.Hyperrectangle(;
    low = SVector(-disc.u_max, -disc.u_max, -disc.u_max),
    high = SVector(disc.u_max, disc.u_max, disc.u_max),
)
@time carving = RP.infeasible_cells(
    X_box,
    disc.state_grid,
    obstacle,
    geometry;
    obstacle_name = "lift constraint",
)
removed = carving.cells
margin = RP.swing_foot_deviation_bound(geometry, MP.get_h(disc.state_grid) ./ 2)
println("Lipschitz carving margin: ", round(100 * margin; digits = 1), " cm")
X = UT.set_minus(X_box, removed)

# How big the abstraction will be, before it is computed: its cost is linear in the number of
# `(cell, input)` pairs, and the number of cells is the product of the cell counts along the five
# axes, so a "slightly" wider state box can multiply it by a lot (see `RP.abstraction_size`).
size = RP.abstraction_size(
    X_box,
    disc.state_grid,
    carving;
    ninputs = RP.input_count(disc; one_joint_per_step),
)
println(
    "State space: ",
    size.cells,
    " cells (",
    size.cells_in_box,
    " before carving), ",
    size.pairs,
    " (cell, input) pairs: about ",
    round(size.transitions / 1e6; digits = 1),
    " M transitions and ",
    round(size.minutes; digits = 1),
    " min on ",
    Threads.nthreads(),
    " thread(s)",
)

# ------------------------------------------------------------
# 3) Abstraction with local growth bounds
# ------------------------------------------------------------

@time lgb = LGA.LocalGrowthBound(f_sym, x_sym, u_sym, X, U_box, disc.tstep)
concrete_system = ST.get_system(lgb)   # discrete-time (RK4) system on the carved domain

# Concrete start of the closed-loop run below: the configuration `q0` with the first coordinate in the
# middle of the initial range, snapped to the grid so that the run starts where the abstraction does
# (an error, naming the constraint that did it, if its cell was carved out).
x0 = RP.start_state(
    disc.state_grid,
    q0,
    X_box;
    velocity = initial_velocity,
    removed = carving,
    geometry,
)
println("x0 = ", x0, "  (ṗ₁ at x0: ", lgb.dynamics.f(x0, SVector(0.0, 0.0, 0.0))[1], ")")
println("swing foot at x0: ", RP.swing_foot_position(SVector{4}(x0[2:5]), geometry))

# `build_abstraction` runs the Dionysos `UniformGridAbstraction` optimizer with the reach sets of
# `lgb`; the returned optimizer has the abstraction cached.
#
# **Computing it on a server.** The abstraction is by far the most expensive part, and it does not
# depend on the control problem. With `ABSTRACTION_FILE` set, its transitions are cached in that file:
# loaded if the file exists (no computation), otherwise computed and written there. With
# `ABSTRACTION_ONLY = true` the run stops once the abstraction is available (nothing else is done). The
# grids, the state box and the other parameters above are stored with the transitions and checked at
# loading: a file computed for other parameters is refused. For a **Slurm array**, the transitions are
# computed in chunks by `EXECUTION_BACKEND = SY.SlurmArrayBackend(nchunks, chunk_id, outdir,
# :contiguous, true)` (one run per chunk), then `MERGE_CHUNKS = (; outdir, nchunks)` reads them
# instead of computing and `ABSTRACTION_FILE` receives the merged result. See
# `problems/BipedRobot/4D_model_underactuated/slurm_array`.
abstraction_file = @isdefined(ABSTRACTION_FILE) ? ABSTRACTION_FILE : nothing
abstraction_only = @isdefined(ABSTRACTION_ONLY) ? ABSTRACTION_ONLY : false
# An outdated `LocalGrowthAbstraction.jl` would pass `abstraction_file` on as an optimizer attribute that
# nothing reads, and silently compute the abstraction: refuse it.
abstraction_file === nothing ||
    isdefined(LGA, :merge_chunks_to_file) ||
    error(
        "$(PROBLEM_DIR)/LocalGrowthAbstraction.jl is too old for ABSTRACTION_FILE (no `merge_chunks_to_file`). " *
        "Dionysos is loaded from $(pathof(Dionysos)): `PROBLEM_DIR` is relative to it.",
    )
execution_backend = @isdefined(EXECUTION_BACKEND) ? EXECUTION_BACKEND : SY.ThreadedBackend()
merge_chunks = @isdefined(MERGE_CHUNKS) ? MERGE_CHUNKS : nothing
# Merging the chunks into `ABSTRACTION_FILE` does not keep the abstraction in memory (it would not fit,
# that is the point): the run has to stop there, and the next one loads the file.
merge_chunks === nothing ||
    abstraction_file === nothing ||
    abstraction_only ||
    error("`MERGE_CHUNKS` writes `ABSTRACTION_FILE` and stops: set `ABSTRACTION_ONLY = true`")
abstraction_info = (;
    state_origin = collect(MP.get_origin(disc.state_grid)),
    state_h = collect(MP.get_h(disc.state_grid)),
    input_h = collect(MP.get_h(disc.input_grid)),
    tstep = disc.tstep,
    x_lb = collect(x_lb),
    x_ub = collect(x_ub),
    ninputs = RP.input_count(disc; one_joint_per_step),
    obstacle = obstacle === nothing ? nothing : length(obstacle),
)
@time optimizer = LGA.build_abstraction(
    lgb;
    state_grid = disc.state_grid,
    input_grid = disc.input_grid,
    print_level = 1,
    use_implicit_mapping = true,
    execution_backend,
    abstraction_file,
    abstraction_info,
    merge_chunks,
    extra_attributes...,
)
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
if abstract_system === nothing
    println("abstraction written to ", abstraction_file, " (not kept in memory)")
else
    println("transitions: ", SY.get_n_transitions(abstract_system))
end
println("local growth bound statistics: ", LGA.stats(lgb))

if !abstraction_only

    # ------------------------------------------------------------
    # 4) One footstep: reach the target foothold without falling, with any velocity
    # ------------------------------------------------------------

    # Reach the foothold from `q0`, with any velocity at the start and at the end (by default): the
    # over-approximation is nondeterministic, so the synthesis works against the worst case (no
    # deterministic-only cost is given).
    step_pb = RP.step_problem(
        concrete_system,
        disc.state_grid;
        q0 = q0,
        foothold = foothold,
        geometry = geometry,
        initial_velocity = initial_velocity,
        target_velocity = target_velocity,
        foothold_tolerance = foothold_tolerance,
    )

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
    @time MOI.optimize!(optimizer)
    println("footstep success: ", MOI.get(optimizer, MOI.RawOptimizerAttribute("success")))

    # The abstraction is a certified over-approximation (up to the boundary convention of
    # `LocalGrowthBound`), so a successful synthesis is guaranteed; the concrete closed loop below
    # is the end-to-end check.
    traj = RP.closed_loop(optimizer, step_pb, x0; ankle_name = "p₁")

    # ------------------------------------------------------------
    # 5) Visualization
    # ------------------------------------------------------------

    # Define `RUN_PLOTS = false` before including this file to skip it. `RP.plot_footstep` also draws
    # nothing when the controller has no action at `x0` already (the default outcome of the experiment
    # without a controller): `traj` has no input then. The figures and the animation go next to this
    # file, or to the directory `OUTPUT_DIR` if it is defined.
    output_dir = @isdefined(OUTPUT_DIR) ? OUTPUT_DIR : @__DIR__
    mkpath(output_dir)
    if @isdefined(RUN_PLOTS) ? RUN_PLOTS : true
        RP.plot_footstep(
            traj,
            disc.state_grid,
            removed;
            x_lb,
            x_ub,
            tstep = disc.tstep,
            geometry,
            obstacle,
            foothold,
            file_prefix = joinpath(output_dir, "biped_4d_normal_form_$(scenario_name)"),
            ankle_label = "p₁ (N m s)",
            dynamics = lgb.dynamics.f,      # re-integrate each step for a smoother animation
            subframes = 10,
            title = "Biped in normal form — footstep ($(scenario_name))",
        )
    end # run_plots
end # !abstraction_only
