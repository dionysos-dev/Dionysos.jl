# Deploy the certified 4-D biped footstep controller on the control server.
#
# The controller is produced by `examples/BipedRobot/export_controller.jl`, which
# saves it under the key "controller" -- what `load_server_controller` reads.
# Nothing of the abstraction is needed here: a compacted controller carries only
# its state mapping and decision table, so this starts in seconds.
#
#   julia --project=control_server control_server/scripts/deploy_biped_footstep.jl
#
# Environment:
#   BIPED_SCENARIO    scenario name, default "step" (matches the example)
#   BIPED_CONTROLLER  "slew" (default) or "plain"
#   BIPED_PORT        TCP port, default 5000
#
# Prefer "slew": its commands respect the one-notch-per-step velocity limit, so a
# real motor controller can track them. "plain" reaches the same foothold in 27
# steps instead of 30, but jumps a full 1.0 rad/s between consecutive commands.
#
# The robot sends the 4 joint angles and receives the 4 joint velocities to
# command, holding them for `tstep = 0.1 s`. Talk to it with
# `biped_footstep_client.jl`, which plays that loop against a simulated plant.

include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer

const CSD = ControlServer.ControlServerDeployment
const SR = ControlServer.ServerRuntime

const SCENARIO = get(ENV, "BIPED_SCENARIO", "step")
const VARIANT = get(ENV, "BIPED_CONTROLLER", "slew")
const PORT = parse(Int, get(ENV, "BIPED_PORT", "5000"))

VARIANT in ("slew", "plain") ||
    error("BIPED_CONTROLLER must be \"slew\" or \"plain\", got \"$VARIANT\"")

suffix = VARIANT == "slew" ? "_controller_slew" : "_controller"
controller_file = joinpath(
    @__DIR__,
    "..",
    "..",
    "examples",
    "BipedRobot",
    "biped_4d_footstep_$(SCENARIO)$(suffix).jld2",
)

isfile(controller_file) || error(
    "no controller at $controller_file -- generate it first with\n" *
    "    julia -t auto --project=test examples/BipedRobot/export_controller.jl",
)

println(
    "loading $(VARIANT) controller for scenario \"$(SCENARIO)\" ",
    "($(round(filesize(controller_file) / 1024 / 1024; digits = 1)) MB)",
)
server_controller = CSD.load_server_controller(controller_file)
println("controller ready, memory = ", server_controller.x)

# 4 joint angles in, 4 joint velocities out.
result = SR.start_control_server(
    server_controller;
    port = PORT,
    log_data = true,
    received_data_size = 4,
    state_to_vector = nothing,
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
