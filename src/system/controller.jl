# Controllers have two orthogonal axes:
#  - discrete vs continuous (the *level* they act on: automaton states vs concrete
#    states) — expressed by the type hierarchy, used to type solver fields;
#  - static vs dynamic (stateless feedback map vs internal memory) — expressed by
#    the `controller_kind` trait, used to dispatch the closed-loop simulation.
abstract type AbstractController end
abstract type AbstractDiscreteController <: AbstractController end
abstract type AbstractContinuousController <: AbstractController end

"""
    ControllerKind

Trait: `StaticKind()` for a stateless feedback map, `DynamicKind()` for a
controller with internal memory (protocol `initial_state`/`update_state`).
"""
abstract type ControllerKind end
struct StaticKind <: ControllerKind end
struct DynamicKind <: ControllerKind end

controller_kind(ctrl::AbstractController) =
    error("implement `controller_kind` for $(typeof(ctrl))")

initial_state(ctrl::AbstractController) =
    error("implement `initial_state` for $(typeof(ctrl))")
update_state(ctrl::AbstractController, controller_state, measurement) =
    error("implement `update_state` for $(typeof(ctrl))")
output_control(ctrl::AbstractController, controller_state, measurement) =
    error("implement `output_control` for $(typeof(ctrl))")

domain(ctrl::AbstractController) = error("implement `domain` for $(typeof(ctrl))")
is_defined(ctrl::AbstractController, controller_state, measurement) = true

# --------------- Discrete Static Controller --------------------------

struct DiscreteStaticController{D, C} <: AbstractDiscreteController
    dom::D
    controller_map::C
    randomize::Bool
end

controller_kind(::DiscreteStaticController) = StaticKind()
domain(ctrl::DiscreteStaticController) = ctrl.dom
initial_state(::DiscreteStaticController) = nothing
update_state(::DiscreteStaticController, x, y) = nothing

function is_defined(ctrl::DiscreteStaticController, x, q::Int)
    q in domain(ctrl) || return false
    us = ctrl.controller_map(q)
    return us !== nothing && !isempty(us)
end

function output_control(ctrl::DiscreteStaticController, x, q::Int)
    is_defined(ctrl, x, q) || return nothing
    us = ctrl.controller_map(q)
    return ctrl.randomize ? rand(us) : first(us)
end

# --------------- Discrete Dynamic Controller --------------------------

struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

struct DiscreteDynamicController{XD, D, G, H} <: AbstractDiscreteController
    x0::XD
    dom::D
    statemap::G
    outputmap::H
    randomize::Bool
end

controller_kind(::DiscreteDynamicController) = DynamicKind()
domain(ctrl::DiscreteDynamicController) = ctrl.dom
initial_state(ctrl::DiscreteDynamicController) = ctrl.x0

function update_state(ctrl::DiscreteDynamicController, x, y)
    return ctrl.statemap(x, y)
end

function is_defined(ctrl::DiscreteDynamicController, x, y)
    (x, y) in domain(ctrl) || return false
    us = ctrl.outputmap(x, y)
    return us !== nothing && !(us isa AbstractVector && isempty(us))
end

function output_control(ctrl::DiscreteDynamicController, x, y)
    is_defined(ctrl, x, y) || return nothing
    us = ctrl.outputmap(x, y)
    return us isa AbstractVector ? (ctrl.randomize ? rand(us) : first(us)) : us
end
