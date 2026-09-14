# Deploy the certified 4-D biped footstep controller on the control server.
#
# The controller is produced by `examples/BipedRobot/export_controller.jl`, which
# saves it under the key "controller" -- what `load_server_controller` reads.
# Nothing of the abstraction comes with it: a compacted controller carries only
# its state mapping and decision table, so this starts in seconds.
#
#   julia --project=control_server control_server/scripts/deploy_biped_footstep.jl
#
# Point `controller_file` at the plain variant to compare: it reaches the same
# foothold in 27 steps instead of 30, but jumps a full 1.0 rad/s between
# consecutive commands, where the slew-limited one stays within 0.5 -- the
# acceleration limit a real motor controller can actually track.
#
# `biped_footstep_client.jl` stands in for the robot; run it against this to
# check a deployment before the real thing.

include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer
import JLD2

const CSD = ControlServer.ControlServerDeployment
const SR = ControlServer.ServerRuntime

const PORT = 5000
const DEGREES = true  # units on the wire. Dionysos is radians throughout
const controller_file = joinpath(@__DIR__, "biped_4d_footstep_step_controller_slew.jld2")

isfile(controller_file) || error(
    "no controller at $(abspath(controller_file)) -- generate it with\n" *
    "    julia -t auto --project=test examples/BipedRobot/export_controller.jl\n" *
    "then copy the .jld2 next to this script.",
)

# Loading is ~166 MB and a robot session re-includes this file often. Keep the
# controller across reloads -- but cache the *Dionysos* controller, not the
# server wrapper: the server mutates the wrapper's `.x` as the session runs and
# takes its initial memory from whatever it holds at startup, so reusing a
# wrapper would silently start the next session from the last one's memory.
if !@isdefined(concrete_controller)
    println(
        "loading ",
        basename(controller_file),
        " (",
        round(filesize(controller_file) / 1024 / 1024; digits = 1),
        " MB)",
    )
    const concrete_controller = JLD2.load(controller_file, "controller")
else
    println("controller already in memory, skipping the reload")
end

server_controller = CSD.to_server_controller(concrete_controller)

# The controller is in radians -- joint angles in, rad/s out -- but LabVIEW
# speaks degrees. Convert at the boundary, leaving the memory to the server:
# `f` and `g` take it as their first argument, so wrapping them does not touch
# it.
if DEGREES
    inner = server_controller
    server_controller = ControlServer.Controller4Server.Controller(
        inner.x,
        (mem, y) -> inner.f(mem, deg2rad.(y)),
        (mem, y) -> begin
            u = inner.g(mem, deg2rad.(y))
            return u === nothing ? nothing : rad2deg.(u)
        end,
    )
end

println(
    "controller ready, memory = ",
    server_controller.x,
    ", wire units = ",
    DEGREES ? "degrees" : "radians",
)

# 4 joint angles in, 4 joint velocities out.
result = SR.start_control_server(
    server_controller;
    port = PORT,
    log_data = true,
    received_data_size = 4,
    state_to_vector = nothing,
    # The controller commands a velocity to hold for one `tstep`. Its slew bound
    # is `du` per call, so it is the intended acceleration limit only if the
    # client keeps that rate -- warn early if it does not.
    expected_dt = 0.1,
)

if result !== nothing
    t, measurements, controls, _ = result
    n = length(t)
    println("session finished: $n packets")
    if n > 1
        println("average dt: ", round((t[n] - t[1]) / (n - 1); digits = 4), " s")
        slew =
            maximum(maximum(abs.(controls[:, k + 1] - controls[:, k])) for k in 1:(n - 1))
        println("max |Δu| along the session: ", round(slew; digits = 3), " rad/s")
    end
end

# One panel per joint: what the robot reported against what it was told to do.
# `Plots` is not a dependency of this environment on purpose -- it resolves from
# the default one through the stacked load path. Guarded so a machine without it
# loses the figure rather than the session that has already finished.
if result !== nothing && n > 0
    if Base.find_package("Plots") === nothing
        @info "Plots is not in the load path; skipping the session figure"
    else
        using Plots
        unit = DEGREES ? "deg" : "rad"
        fig = plot(; layout = (2, 2), size = (1200, 900))
        for i in 1:4
            plot!(fig[i], t[1:n], measurements[i, 1:n]; label = "θ$i ($unit)")
            plot!(fig[i], t[1:n], controls[i, 1:n]; label = "u$i ($unit/s)")
            xlabel!(fig[i], "time (s)")
        end
        display(fig)
    end
end
