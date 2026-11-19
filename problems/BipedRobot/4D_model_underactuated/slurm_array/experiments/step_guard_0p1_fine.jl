# The step of `step_guard_0p1.jl` (`tstep = 1/60 s`) on a finer grid along the ankle only: `q₁` in cells of
# 0.05 rad instead of 0.1, and the momentum `p₁` in cells of 0.015 instead of 0.04; the motor angles,
# the box, the obstacle, the inputs and the target are unchanged. At 0.1/0.04 the ankle moves less than
# one cell per step (`|Δp₁| ≤ 0.0185`; `q₁` too, for every input, in half of the cells), so the
# abstraction can keep the ankle in its cell at every step and the worst case does: the domain of
# that run never left the target's range of `q₁`. At `1/20 s` (`step_guard_0p1_t20.jl`) the ankle
# moves on, but the images are wider and the domain halves. Here the cells shrink to the ankle's
# motion instead of the motion growing to the cells.
#
# The parameter file of the Slurm pipeline (`DIONYSOS_PARAMS`, README §2): included before
# `examples/BipedRobot/biped_4d_normal_form.jl`, it sets the variables the example reads. Pass the same
# file to every job (the chunks, the merge and the solve job): the obstacle, the grid and the box
# enter the abstraction.
#
# * Grid: 0.1 rad on the motor angles `q₂, q₃, q₄`, 0.05 rad on the ankle `q₁`, 0.015 N m s on the ankle
#   momentum `p₁`, `tstep = 1/60 s` (the motors move by one cell per step at `du = 6 rad/s`), and five
#   speeds per motor joint, `{0, ±du, ±2du}`: 125 inputs.
# * State box: as in `step_guard_0p1.jl`, `p₁ ∈ [-0.6, -0.08]` and the knee box of size 1.5:
#   11 402 104 cells, 7 120 416 after carving (5.7 times the 0.1/0.04 grid), 890 M `(cell, input)`
#   pairs (about 2.5 ms each on a Rome core, measured on the earlier runs: about 620 core-hours).
#   On the example's smaller box (270 232 cells at this grid), the finer ankle has 1.36 times as many
#   transitions per cell as 0.1/0.04, hence about 12 billion here (1 285 per cell at 0.1/0.04), so a
#   solve of about 600 GB at 50 bytes per transition: a 1 TB node of `keira`. On that smaller box, the ankle keeps
#   both its cells for 38 % of the `(cell, input)` pairs instead of 85 %, and every input lets it stay
#   in 5 % of the cells instead of 49 %.
# * Ground, target and obstacle: as in `step_guard_0p1.jl`. The target has 17 646 cells (519 angle
#   cells), none of them carved by the obstacle; the start `q0` is kept (`p₁ = -0.405`, the cell of
#   width 0.015 that contains `-0.4`).
#
# The start, the target momentum and the tolerance do not enter the abstraction: change them in a
# `DIONYSOS_SOLVE_PARAMS` file for the solve job.

using StaticArrays
import LazySets

DISC_KWARGS = (; dx = 0.1, dq1 = 0.05, dq̇1 = 0.015, tstep = 1 / 60, input_steps = 2)
X_BOUNDS = (SVector(-0.6, -1.5, 0.0, -1.5, -1.5), SVector(-0.08, 1.5, 1.5, 1.5, 0.0))
OBSTACLE =
    [LazySets.Hyperrectangle(; low = SVector(-0.135, -0.15), high = SVector(0.139, 0.0))]
GROUND_AMBIGUOUS = :keep
FOOTHOLD_TOLERANCE = 0.0
