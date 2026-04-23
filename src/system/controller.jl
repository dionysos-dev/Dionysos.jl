abstract type AbstractController end
abstract type AbstractDiscreteController <: AbstractController end
abstract type AbstractContinuousController <: AbstractController end

initial_state(ctrl::AbstractController) = error("Not implemented")
update_state(ctrl::AbstractController, x, y) = error("Not implemented")
output_control(ctrl::AbstractController, x, y) = error("Not implemented")

domain(ctrl::AbstractController) = error("Not implemented")
is_defined(ctrl::AbstractController, x, y) = true

state_domain(ctrl::AbstractController) = nothing
input_domain(ctrl::AbstractController) = nothing

# --------------- Discrete Static Controller --------------------------

struct DiscreteStaticController{D, C} <: AbstractDiscreteController
    dom::D
    controller_map::C
    randomize::Bool
end

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
