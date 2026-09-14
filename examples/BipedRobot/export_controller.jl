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
# Usage — the scenario is selected exactly as in the example itself, and the
# state-grid / input-grid sizes make the run take several minutes:
#
# ```
# julia -t auto --project=test examples/BipedRobot/export_controller.jl
# BIPED_SCENARIO=wall julia -t auto --project=test examples/BipedRobot/export_controller.jl
# ```
#
# `-t auto` matters: the example asks for `SY.ThreadedBackend()`, but that falls
# back to a sequential build when Julia runs single-threaded, which roughly
# doubles the abstraction time.
#
# The `.jld2` files are deliberately untracked (`*.jld2` is in `.gitignore`) —
# they are regenerable artifacts, and committing them is what bloated this
# repository's history in the first place.

import JLD2

# Runs the whole example, including its `.png` / `.gif` output.
include(joinpath(@__DIR__, "biped_4d_velocity.jl"))

# `controller`, `slew_controller`, `scenario_name`, `disc`, `x0` and
# `discrete_time_system` all come from the example.
plain_file = joinpath(@__DIR__, "biped_4d_footstep_$(scenario_name)_controller.jld2")
slew_file = joinpath(@__DIR__, "biped_4d_footstep_$(scenario_name)_controller_slew.jld2")

# `disc` is a NamedTuple of (state_grid, input_grid, tstep, u_max, du).
metadata = (;
    scenario = String(scenario_name),
    tstep = disc.tstep,
    du = disc.du,
    u_max = disc.u_max,
    source = "examples/BipedRobot/biped_4d_velocity.jl",
)

# `compact_controller` drops the abstraction's transition relation, which a
# finished controller never reads again -- ~39 M transitions here, and the
# difference between a file you can ship and one you cannot.
JLD2.jldsave(
    plain_file;
    controller = SY.compact_controller(controller),
    metadata = metadata,
    variant = "plain",
)
JLD2.jldsave(
    slew_file;
    controller = SY.compact_controller(slew_controller),
    metadata = metadata,
    variant = "slew_rate_limited",
)

for file in (plain_file, slew_file)
    println("saved ", file, " (", round(filesize(file) / 1024 / 1024; digits = 2), " MB)")
end

# Round-trip the files: one that cannot be read back is not a controller. The
# closed loop is the real consumer, so replay it rather than poking the
# controller protocol directly.
for (file, reference) in ((plain_file, controller), (slew_file, slew_controller))
    loaded = JLD2.load(file, "controller")
    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        loaded,
        x0,
        400;
        stopping = reached,
    )
    xs = collect(ST.states(traj))
    println(
        basename(file),
        ": type round-trips = ",
        typeof(loaded) === typeof(reference),
        ", replayed ",
        length(xs) - 1,
        " steps, reached = ",
        reached(xs[end]),
    )
end
