# ============================================================
# Controller protocol.
#
# Controllers have two orthogonal axes:
#  - discrete vs continuous (the *level* they act on: automaton states vs concrete
#    states) — expressed by the type hierarchy, used to type solver fields;
#  - static vs dynamic (stateless feedback map vs internal memory) — expressed by
#    the `controller_kind` trait, used to dispatch the closed-loop simulation.
#
# Protocol (argument order is always `(controller, memory, measurement)`):
#  - `controller_kind(ctrl) -> StaticKind() | DynamicKind()`  (required)
#  - `output_control(ctrl, mem, y) -> u` or `nothing` when undefined  (required)
#  - `initial_state(ctrl) -> mem₀` and `update_state(ctrl, mem, y) -> mem′`
#    (required for dynamic controllers; static controllers inherit `nothing`)
#  - `is_defined(ctrl, mem, y) -> Bool` (defaults to `true`; the authoritative
#    membership test — `domain` is optional introspection and may be `nothing`)
#
# Serializability convention: a controller a user may want to save (JLD2) must be
# plain data — tables, sets, numbers, or struct callables (see
# `src/utils/functions.jl`) — never closures.
# ============================================================

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

# Trait-based defaults: static controllers are memoryless, so only dynamic
# controllers must implement the memory half of the protocol.
initial_state(ctrl::AbstractController) = _initial_state(controller_kind(ctrl), ctrl)
_initial_state(::StaticKind, ctrl) = nothing
_initial_state(::DynamicKind, ctrl) = error("implement `initial_state` for $(typeof(ctrl))")

update_state(ctrl::AbstractController, mem, y) =
    _update_state(controller_kind(ctrl), ctrl, mem, y)
_update_state(::StaticKind, ctrl, mem, y) = nothing
_update_state(::DynamicKind, ctrl, mem, y) =
    error("implement `update_state` for $(typeof(ctrl))")

output_control(ctrl::AbstractController, mem, y) =
    error("implement `output_control` for $(typeof(ctrl))")

# `nothing` means "no domain information"; `is_defined` is the authoritative test.
domain(ctrl::AbstractController) = nothing
is_defined(ctrl::AbstractController, mem, y) = true

# --------------- Control table --------------------------

"""
    ControlTable(nstates::Int)

Plain-data map from an abstract state `q ∈ 1:nstates` to the list of admissible
abstract inputs, the canonical `controller_map` of a
[`DiscreteStaticController`](@ref). Callable: `table(q)` returns the input list;
fill with `add_control!` / `set_control!`.
"""
mutable struct ControlTable
    U::Vector{Vector{Int}}
end

ControlTable(nstates::Int) = ControlTable([Int[] for _ in 1:nstates])

(C::ControlTable)(q::Int) = C.U[q]

"Append input `u` to the admissible inputs of state `q`."
add_control!(C::ControlTable, q::Int, u::Int) = push!(C.U[q], u)

"Make `u` the unique admissible input of state `q`."
function set_control!(C::ControlTable, q::Int, u::Int)
    empty!(C.U[q])
    push!(C.U[q], u)
    return u
end

# --------------- Discrete Static Controller --------------------------

"""
    DiscreteStaticController(dom, controller_map, randomize)

Static feedback on abstract states: `controller_map(q)` returns the admissible
inputs of `q` (e.g. a [`ControlTable`](@ref)), `dom` the set of controlled states.
"""
struct DiscreteStaticController{D, C} <: AbstractDiscreteController
    dom::D
    controller_map::C
    randomize::Bool
end

controller_kind(::DiscreteStaticController) = StaticKind()
domain(ctrl::DiscreteStaticController) = ctrl.dom

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

"""
    DiscreteDynamicController(x0, dom, statemap, outputmap, randomize)

Dynamic feedback on abstract states with memory `x` starting at `x0`:
`statemap(x, y)` updates the memory, `outputmap(x, y)` returns the control (or
the list of admissible controls), `dom` contains the valid `(x, y)` pairs.
Closure-backed and therefore not serializable — prefer
[`AutomatonMemoryController`](@ref) for controllers meant to be saved.
"""
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

# --------------- Automaton-Memory Controller --------------------------

"""
    AutomatonMemoryController

Finite-memory controller over abstract states, fully table-backed — the
serializable alternative to a closure-based [`DiscreteDynamicController`](@ref)
for specification-automaton memory (e.g. co-safe LTL):

- the memory follows the specification automaton:
  `qa′ = step_map[(qa, label)]` with `label = label_of_state[qs]`
  (`default_label` for states outside the labeled range);
- the control comes from a controller synthesized on the product automaton,
  looked up through the product-state table `pid[(qs, qa)]`.
"""
struct AutomatonMemoryController{C <: AbstractController} <: AbstractDiscreteController
    x0::Int
    label_of_state::Vector{Int}
    default_label::Int
    step_map::Dict{Tuple{Int, Int}, Int}
    pid::Dict{Tuple{Int, Int}, Int}
    product_controller::C
end

controller_kind(::AutomatonMemoryController) = DynamicKind()
initial_state(ctrl::AutomatonMemoryController) = ctrl.x0

_label_id(ctrl::AutomatonMemoryController, qs::Int) =
    1 <= qs <= length(ctrl.label_of_state) ? ctrl.label_of_state[qs] : ctrl.default_label

update_state(ctrl::AutomatonMemoryController, qa, qs::Int) =
    get(ctrl.step_map, (qa, _label_id(ctrl, qs)), nothing)

function is_defined(ctrl::AutomatonMemoryController, qa, qs::Int)
    p = get(ctrl.pid, (qs, qa), nothing)
    p === nothing && return false
    return is_defined(ctrl.product_controller, nothing, p)
end

function output_control(ctrl::AutomatonMemoryController, qa, qs::Int)
    p = get(ctrl.pid, (qs, qa), nothing)
    p === nothing && return nothing
    return output_control(ctrl.product_controller, nothing, p)
end

# --------------- Affine Controller --------------------------

"""
    AffineController(map::MathematicalSystems.AffineMap)

Static state feedback `u = A·x + c` (plain data). Wrap a transition-synthesis
result with `as_controller` to simulate it in a closed loop.
"""
struct AffineController{M <: MS.AffineMap} <: AbstractContinuousController
    map::M
end

controller_kind(::AffineController) = StaticKind()
output_control(ctrl::AffineController, mem, x) = MS.apply(ctrl.map, x)

# --------------- Funnel Controller --------------------------

"""
    FunnelController(kappas, ellipsoids)

Time-indexed funnel-tracking controller, the simulatable form of an ellipsoidal
certification chain `(E_1, κ_1, …, E_K, κ_K, E_{K+1})`: the memory is the step
index `k`, the output at step `k` is `κ_k(x)` (an absolute-coordinates
`MathematicalSystems.AffineMap`), and the controller is defined while the measured
state lies in the certified ellipsoid `E_k`. The certificate guarantees that a
defined step stays defined: `x ∈ E_k ⟹ x⁺ ∈ E_{k+1}`. Plain data — serializable.
"""
struct FunnelController{TK, TE} <: AbstractContinuousController
    kappas::Vector{TK}
    ellipsoids::Vector{TE}

    function FunnelController(kappas::Vector{TK}, ellipsoids::Vector{TE}) where {TK, TE}
        length(ellipsoids) == length(kappas) + 1 ||
            error("need length(ellipsoids) == length(kappas) + 1 (E_1..E_{K+1}, κ_1..κ_K)")
        return new{TK, TE}(kappas, ellipsoids)
    end
end

controller_kind(::FunnelController) = DynamicKind()
initial_state(::FunnelController) = 1
update_state(ctrl::FunnelController, k, x) = k + 1
domain(ctrl::FunnelController) = ctrl.ellipsoids

function is_defined(ctrl::FunnelController, k, x)
    1 <= k <= length(ctrl.kappas) || return false
    return collect(x) ∈ ctrl.ellipsoids[k]
end

function output_control(ctrl::FunnelController, k, x)
    is_defined(ctrl, k, x) || return nothing
    return MS.apply(ctrl.kappas[k], x)
end
