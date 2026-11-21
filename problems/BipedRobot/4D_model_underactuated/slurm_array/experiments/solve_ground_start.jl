# Solve parameters (`DIONYSOS_SOLVE_PARAMS` of `solve.jl`, not an abstraction): the step from every
# configuration with the swing foot on the ground around `(-0.2, 0)`, `x` every centimetre from -0.25 to
# -0.15, with any momentum, instead of the single configuration of `x0` (whose swing foot is 2.8 cm in
# the air). For the abstraction of `step_guard_0p1_fine.jl`. Its saved domain (which does not depend on
# the start: the synthesis runs to its fixed point) already has 762 of these 95,064 start cells, 230
# with the foot exactly at `(-0.2, 0)`, with momenta from -0.465 to -0.24.
START_FOOTHOLDS = [SVector(x, 0.0) for x in -0.25:0.01:-0.15]
