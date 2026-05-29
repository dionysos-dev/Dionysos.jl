using StaticArrays
using Dionysos

const ST = Dionysos.System

struct LiftedRobotController{C} <: ST.AbstractContinuousController
    controller6::C
end

ST.initial_state(ctrl::LiftedRobotController) =
    ST.initial_state(ctrl.controller6)

# 8D measured state -> 6D reduced state
@inline function project_state_8_to_6(x8)
    return SVector{6,Float64}(
        x8[1], # q3
        x8[2], # q4
        x8[3], # q5
        x8[5], # dq3, adapt index if velocity layout differs
        x8[6], # dq4
        x8[7], # dq5
    )
end

# 3D symbolic input -> 4D real input
@inline function lift_input_3_to_4(u3)
    return SVector{4,Float64}(
        u3[1],
        u3[2],
        u3[3],
        0.0, # fourth input handled by low-level knee controller
    )
end

function ST.is_defined(ctrl::LiftedRobotController, q, y8)
    y6 = project_state_8_to_6(y8)
    return ST.is_defined(ctrl.controller6, q, y6)
end

function ST.output_control(ctrl::LiftedRobotController, q, y8)
    y6 = project_state_8_to_6(y8)
    u3 = ST.output_control(ctrl.controller6, q, y6)
    u3 === nothing && return nothing
    return lift_input_3_to_4(u3)
end

function ST.update_state(ctrl::LiftedRobotController, q, y8)
    y6 = project_state_8_to_6(y8)
    return ST.update_state(ctrl.controller6, q, y6)
end