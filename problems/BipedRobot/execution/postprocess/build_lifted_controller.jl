using JLD2
using Dionysos
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "6D_model", "two_step_walking_controller.jl"))
include(joinpath(@__DIR__, "..", "..", "6D_model", "lifted_robot_controller.jl"))

const CONTROLLER_6D_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_two_step_controller_CM.jld2")

const LIFTED_CONTROLLER_FILE =
    joinpath(@__DIR__, "..", "..", "out", "6D", "robot_two_step_controller_8D4U_CM.jld2")

function load_controller_data(filename::AbstractString)
    return JLD2.jldopen(filename, "r") do file
        return (
            controller = file["controller"],
            tstep = file["tstep"],
        )
    end
end

function save_lifted_controller(filename::AbstractString; controller, tstep)
    mkpath(dirname(filename))

    JLD2.jldopen(filename, "w") do file
        file["controller"] = controller
        file["tstep"] = tstep
        file["state_dimension"] = 8
        file["input_dimension"] = 4
        file["measurement_convention"] =
            "[q1,...,q8,dq1,...,dq8], using q3,q4,q5,dq3,dq4,dq5"
    end

    return filename
end

@info "Loading 6D controller" CONTROLLER_6D_FILE
data = load_controller_data(CONTROLLER_6D_FILE)

lifted_controller = LiftedRobotController(data.controller)

save_lifted_controller(
    LIFTED_CONTROLLER_FILE;
    controller = lifted_controller,
    tstep = data.tstep,
)

@info "Saved lifted 8D/4-input controller" LIFTED_CONTROLLER_FILE