# # Export the certified footstep controllers to `.jld2`
#
# Runs `biped_4d_velocity.jl` and serializes the two controllers it synthesizes
# so they can be deployed without re-running the abstraction:
#
#   * `<scenario>_controller.jld2`      — the plain footstep controller
#   * `<scenario>_controller_slew.jld2` — the slew-rate-limited one
#     (`OPDS.BoundedInputVariation`), the profile a real motor controller can
#     track, and the one to deploy
#
# Each file stores the controller under the key `"controller"`, which is what
# `control_server/src/ControlServerDeployment.jl` reads via
# `JLD2.load(filename, "controller")`.
#
# Run it with:
#
# ```
# julia -t auto --project=test examples/BipedRobot/export_controller.jl
# ```
#
# `-t auto` is not optional in practice: the example asks for
# `SY.ThreadedBackend()`, which silently falls back to a sequential build on one
# thread — 116 s against 28 s for the abstraction here. Expect several minutes
# either way; the grid is ≈ 2.5 M cells.
#
# To export a different scenario, change `SCENARIO` below — `step`, `riccardo`,
# `wall`, `limbo` or `window`, described at the top of `biped_4d_velocity.jl`.
#
# The `.jld2` files are deliberately untracked (`*.jld2` is in `.gitignore`) —
# they are regenerable artifacts, and committing them is what bloated this
# repository's history in the first place. Copy them next to
# `control_server/scripts/deploy_biped_footstep.jl` to deploy.

import JLD2

# ---------------------------------------------------------------------------
SCENARIO = :step   # step, riccardo, wall, limbo, window
SLEW = true        # true: |Δu| ≤ 0.5, 30 steps, the one to deploy
# ---------------------------------------------------------------------------

RUN_SLEW = SLEW    # read by the example: skip that pass (≈ 165 s) when unwanted
RUN_PLOTS = false  # its .png / .gif are the example's job, not this script's

include(joinpath(@__DIR__, "biped_4d_velocity.jl"))

suffix = SLEW ? "_controller_slew" : "_controller"
file = joinpath(@__DIR__, "biped_4d_footstep_$(SCENARIO)$(suffix).jld2")

# `compact_controller` drops the abstraction's transition relation, which a
# finished controller never reads again -- ~39 M transitions here, and the
# difference between a file you can ship and one you cannot.
JLD2.jldsave(
    file;
    controller = SY.compact_controller(SLEW ? slew_controller : controller),
    variant = SLEW ? "slew_rate_limited" : "plain",
    metadata = (;
        scenario = String(SCENARIO),
        tstep = disc.tstep,
        du = disc.du,
        u_max = disc.u_max,
        source = "examples/BipedRobot/biped_4d_velocity.jl",
    ),
)
println("saved ", file, " (", round(filesize(file) / 1024 / 1024; digits = 2), " MB)")

# A file that cannot be read back is not a controller, and the closed loop is its
# real consumer -- so replay it rather than poking the protocol directly.
loaded = JLD2.load(file, "controller")
replay =
    ST.get_closed_loop_trajectory(discrete_time_system, loaded, x0, 400; stopping = reached)
replay_xs = collect(ST.states(replay))
println("replayed ", length(replay_xs) - 1, " steps, reached = ", reached(replay_xs[end]))
