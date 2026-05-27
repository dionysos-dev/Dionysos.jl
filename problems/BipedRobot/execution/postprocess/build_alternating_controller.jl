using Dionysos
using JLD2
using StaticArrays

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

const NOMINAL_CONTROLLER_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_controller_2.jld2")

const ALTERNATING_CONTROLLER_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_alternating_controller_2.jld2")

include(joinpath(@__DIR__, "..", "..", "6D_model", "alternating_walking_controller.jl"))

function load_controller(filename)
    return JLD2.jldopen(filename, "r") do file
        return (
            controller = file["controller"],
            concrete_system = file["concrete_system"],
            control_problem = file["control_problem"],
            tstep = file["tstep"],
        )
    end
end

function save_alternating_controller(
    filename;
    controller,
    concrete_system,
    control_problem,
    tstep,
)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        file["controller"] = controller
        file["concrete_system"] = concrete_system
        file["control_problem"] = control_problem
        return file["tstep"] = tstep
    end
    return filename
end

data = load_controller(NOMINAL_CONTROLLER_FILE)

alternating_controller = build_alternating_walking_controller(data.controller)

save_alternating_controller(
    ALTERNATING_CONTROLLER_FILE;
    controller = alternating_controller,
    concrete_system = data.concrete_system,
    control_problem = data.control_problem,
    tstep = data.tstep,
)

@info "Saved alternating walking controller" ALTERNATING_CONTROLLER_FILE

# 1. compute abstraction
# 2. synthesize nominal one-step controller
# 3. wrap nominal controller into alternating controller
# 4. simulate / deploy alternating controller
