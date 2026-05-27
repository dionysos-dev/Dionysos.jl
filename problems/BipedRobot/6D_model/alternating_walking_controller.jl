using StaticArrays
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

# ============================================================
# Symmetry maps
# Case 2: same target set, stance-relative coordinates
# ============================================================

@inline function Sx_swap(x::SVector{6, <:Real})
    return SVector{6, Float64}(x[2], x[1], x[3], x[5], x[4], x[6])
end

@inline function Su_swap(u::SVector{3, <:Real})
    return SVector{3, Float64}(u[2], u[1], u[3])
end

@inline function Sx_signed(x::SVector{6, <:Real})
    return SVector{6, Float64}(-x[2], -x[1], x[3], -x[5], -x[4], x[6])
end

@inline function Su_signed(u::SVector{3, <:Real})
    return SVector{3, Float64}(-u[2], -u[1], u[3])
end

# Choose one convention here.
const Sx = Sx_swap
const Su = Su_swap

# ============================================================
# Controller state
# ============================================================

Base.@kwdef struct AlternatingControllerState{Y}
    phase::Symbol = :nominal
    inner_state::Y = nothing
    was_in_target::Bool = false
end

# ============================================================
# Alternating controller wrapper
# ============================================================

struct AlternatingWalkingController{C, T} <: ST.AbstractContinuousController
    controller::C
    target::T
end

function ST.initial_state(ctrl::AlternatingWalkingController)
    return AlternatingControllerState(;
        phase = :nominal,
        inner_state = ST.initial_state(ctrl.controller),
        was_in_target = false,
    )
end

ST.domain(ctrl::AlternatingWalkingController) = ST.domain(ctrl.controller)

ST.state_domain(ctrl::AlternatingWalkingController) = ST.state_domain(ctrl.controller)

ST.input_domain(ctrl::AlternatingWalkingController) = ST.input_domain(ctrl.controller)

# ============================================================
# Active measurement/input according to current phase
# ============================================================

@inline function active_y(
    ctrl::AlternatingWalkingController,
    q::AlternatingControllerState,
    y,
)
    return q.phase == :nominal ? y : Sx(y)
end

@inline function active_u(q::AlternatingControllerState, u)
    u === nothing && return nothing
    return q.phase == :nominal ? u : Su(u)
end

# ============================================================
# Controller interface
# Convention:
#   q = controller internal state
#   y = measured plant state
# ============================================================

function ST.is_defined(ctrl::AlternatingWalkingController, q::AlternatingControllerState, y)
    ya = active_y(ctrl, q, y)
    return ST.is_defined(ctrl.controller, q.inner_state, ya)
end

function ST.output_control(
    ctrl::AlternatingWalkingController,
    q::AlternatingControllerState,
    y,
)
    ya = active_y(ctrl, q, y)

    ua = ST.output_control(ctrl.controller, q.inner_state, ya)

    return active_u(q, ua)
end

function ST.update_state(
    ctrl::AlternatingWalkingController,
    q::AlternatingControllerState,
    y,
)
    in_target = y ∈ ctrl.target

    new_phase =
        in_target && !q.was_in_target ? (q.phase == :nominal ? :converse : :nominal) :
        q.phase

    ya = new_phase == :nominal ? y : Sx(y)

    new_inner_state = ST.update_state(ctrl.controller, q.inner_state, ya)

    return AlternatingControllerState(;
        phase = new_phase,
        inner_state = new_inner_state,
        was_in_target = in_target,
    )
end

# ============================================================
# Builder
# ============================================================

function build_alternating_walking_controller(controller)
    target_low = SVector{6, Float64}(-12π / 180, 7π / 180, 8π / 180, -0.75, -0.30, -0.30)

    target_high = SVector{6, Float64}(-8π / 180, 9π / 180, 12π / 180, 0.30, 0.75, 0.75)

    target = UT.HyperRectangle(target_low, target_high)

    margin = SVector{6, Float64}(1π / 180, 1π / 180, 1π / 180, 0.05, 0.05, 0.05)

    target = UT.HyperRectangle(target.lb .- margin, target.ub .+ margin)

    return AlternatingWalkingController(controller, target)
end
