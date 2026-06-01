using StaticArrays
using JLD2
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

include(joinpath(@__DIR__, "..", "..", "6D_model", "two_step_walking_controller_CM.jl"))

const CONTROLLER_1_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_one_step_controller_1_CM.jld2")

const CONTROLLER_2_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_one_step_controller_2_CM.jld2")

const TWO_STEP_CONTROLLER_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_two_step_controller_CM.jld2")

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

function save_controller(filename; controller, concrete_system, control_problems, tstep)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        file["controller"] = controller
        file["concrete_system"] = concrete_system
        file["control_problems"] = control_problems
        file["control_problem"] = control_problems[1]
        return file["tstep"] = tstep
    end
    return filename
end

data1 = load_controller(CONTROLLER_1_FILE)
data2 = load_controller(CONTROLLER_2_FILE)

two_step_controller = TwoStepWalkingController(
    data1.controller,
    data2.controller,
    data1.control_problem.target_set,
    data2.control_problem.target_set,
)

save_controller(
    TWO_STEP_CONTROLLER_FILE;
    controller = two_step_controller,
    concrete_system = data1.concrete_system,
    control_problems = (data1.control_problem, data2.control_problem),
    tstep = data1.tstep,
)

@info "Saved two-step walking controller" TWO_STEP_CONTROLLER_FILE
