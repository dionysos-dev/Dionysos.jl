using Dionysos

const ST = Dionysos.System

Base.@kwdef struct TwoStepControllerState{Q1, Q2}
    phase::Int = 1
    q1::Q1 = nothing
    q2::Q2 = nothing
    was_in_target::Bool = false
end

struct TwoStepWalkingController{C1, C2, T1, T2} <: ST.AbstractContinuousController
    controller1::C1
    controller2::C2
    target1::T1
    target2::T2
end

ST.controller_kind(::TwoStepWalkingController) = ST.DynamicKind()

function ST.initial_state(ctrl::TwoStepWalkingController)
    return TwoStepControllerState(;
        phase = 1,
        q1 = ST.initial_state(ctrl.controller1),
        q2 = ST.initial_state(ctrl.controller2),
        was_in_target = false,
    )
end

function active_controller(ctrl::TwoStepWalkingController, q::TwoStepControllerState)
    return q.phase == 1 ? ctrl.controller1 : ctrl.controller2
end

function active_state(q::TwoStepControllerState)
    return q.phase == 1 ? q.q1 : q.q2
end

function active_target(ctrl::TwoStepWalkingController, q::TwoStepControllerState)
    return q.phase == 1 ? ctrl.target1 : ctrl.target2
end

function switch_phase(phase::Int)
    return phase == 1 ? 2 : 1
end

function ST.is_defined(ctrl::TwoStepWalkingController, q::TwoStepControllerState, y)
    c = active_controller(ctrl, q)
    qc = active_state(q)
    return ST.is_defined(c, qc, y)
end

function ST.output_control(ctrl::TwoStepWalkingController, q::TwoStepControllerState, y)
    c = active_controller(ctrl, q)
    qc = active_state(q)
    return ST.output_control(c, qc, y)
end

function ST.update_state(ctrl::TwoStepWalkingController, q::TwoStepControllerState, y)
    in_target = y ∈ active_target(ctrl, q)

    can_switch = in_target && !q.was_in_target

    new_phase = can_switch ? switch_phase(q.phase) : q.phase

    q1_new = q.q1
    q2_new = q.q2

    if q.phase == 1
        q1_new = ST.update_state(ctrl.controller1, q.q1, y)
    else
        q2_new = ST.update_state(ctrl.controller2, q.q2, y)
    end

    if can_switch
        @info "Switching walking controller" old_phase = q.phase new_phase = new_phase state =
            y
    end

    return TwoStepControllerState(;
        phase = new_phase,
        q1 = q1_new,
        q2 = q2_new,
        was_in_target = in_target,
    )
end
