# A forward step from x ≈ -0.2 m to the foothold at x = +0.2 m, on the knee box of size 1.5, at 0.1 rad,
# with the cells across the ground kept and the ground between the footholds made solid.
#
# The parameter file of the Slurm pipeline (`DIONYSOS_PARAMS`, README §2): included before
# `examples/BipedRobot/biped_4d_normal_form.jl`, it sets the variables the example reads. Pass the same
# file to every job (the chunks, the merge and the solve job): the obstacle, the grid and the box
# enter the abstraction.
#
# * Grid: 0.1 rad on the four angles, 0.04 N m s on the ankle momentum `p₁`, `tstep = 1/60 s` (the
#   motors move by one cell per step at `du = 6 rad/s`), and five speeds per motor joint,
#   `{0, ±du, ±2du}`: 125 inputs.
# * State box: `p₁ ∈ [-0.6, -0.08]` (a forward step, as in the example) and the knee box of size 1.5
#   for the angles (`q₂ ≥ 0`, `q₄ ≤ 0`): 1 978 032 cells, 1 245 168 after carving, 156 M
#   `(cell, input)` pairs (about 0.5 ms each on one core; about 15 transitions per pair measured on
#   three chunks, hence about 2.3 billion transitions).
# * Ground: `ground_ambiguous = :keep` (the default): the cells whose foot box crosses `z = 0` are
#   kept, since removing them would remove the footholds too.
# * Target: the swing foot exactly on the foothold `(0.2, 0)` (`FOOTHOLD_TOLERANCE = 0`): the cells of
#   `RP.target_cells`, which samples `(q₁, q₂)` 3 times per cell.
# * Obstacle: the rectangle `x ∈ [-0.135, 0.139]`, `z ∈ [-0.15, 0]`. A kept cell whose foot box
#   overlaps it along `x` and reaches `z ≤ 0` is removed, so over that span the swing foot is above the
#   ground in the whole cell (certified, at the sampling times; not between two samples). Its top is
#   the ground: it forbids negative heights there and nothing more. It is the widest such rectangle
#   that keeps every cell of the target set (291 angle cells, whose foot boxes reach down to
#   `x = 0.13903`) and of its mirror image at `(-0.2, 0)` (the starts with the foot on that point: up to
#   `x = -0.13567`; the start `q0` of the example is kept too, its cell is above the ground). Its
#   bottom is below the deepest kept cell (`z_min = -12.85 cm`). Between the rectangle and the
#   footholds the foot can still go below the ground, down to about -9 cm (the profile of
#   `tests/step_guard_sizing.jl`).
#
# The start, the target momentum and the tolerance do not enter the abstraction: change them in a
# `DIONYSOS_SOLVE_PARAMS` file for the solve job. A tolerance, or a denser sampling of the target
# (`nsub = 6` adds cells reaching 4 mm further in), makes the obstacle remove part of the target.

using StaticArrays
import LazySets

DISC_KWARGS = (; dx = 0.1, dq1 = 0.1, dq̇1 = 0.04, tstep = 1 / 60, input_steps = 2)
X_BOUNDS = (SVector(-0.6, -1.5, 0.0, -1.5, -1.5), SVector(-0.08, 1.5, 1.5, 1.5, 0.0))
OBSTACLE =
    [LazySets.Hyperrectangle(; low = SVector(-0.135, -0.15), high = SVector(0.139, 0.0))]
GROUND_AMBIGUOUS = :keep
FOOTHOLD_TOLERANCE = 0.0
