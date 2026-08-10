# Composable trajectory costs (plan.md §3.2-M3). Each term is evaluated *online*
# during a rollout through a functional-accumulator protocol — no trajectory object,
# no per-sample allocation:
#
#     acc = cost_init(term)                 # term-specific accumulator (immutable)
#     acc = cost_step(term, acc, x, u, k)   # once per applied input
#     J   = cost_final(term, acc, xT)       # terminal state closes the account
#
# `CompositeCost` sums its terms. A hand-written `(problem, traj) -> cost` closure
# remains the escape hatch on the generators, but the terms are the recommended path:
# they hoist references once and keep the sampling loop type-stable.

import LazySets
import LinearAlgebra as LA
import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System

export CompositeCost,
    TrackingCost,
    TerminalEllipsoidCost,
    InputEffortCost,
    InputSmoothnessCost,
    DomainPenaltyCost,
    ReachObjectiveCost

abstract type AbstractCostTerm end

cost_init(::AbstractCostTerm) = 0.0
cost_step(::AbstractCostTerm, acc, x, u, k) = acc
cost_final(::AbstractCostTerm, acc, xT) = acc

"""
    CompositeCost(terms...)

Sum of [`AbstractCostTerm`](@ref)s, evaluated online during rollouts.
"""
struct CompositeCost{T <: Tuple}
    terms::T
end
CompositeCost(terms...) = CompositeCost(terms)

"""
    TrackingCost(reference, weights)

Stage cost `Σₖ Σᵢ weights[i]·(x[i] − reference[k][i])²` against a reference state
sequence (hoisted once — never re-collected per sample). Steps beyond the reference
length track its last state.
"""
struct TrackingCost{R, W} <: AbstractCostTerm
    reference::R
    weights::W
end

function cost_step(t::TrackingCost, acc, x, u, k)
    r = t.reference[min(k, length(t.reference))]
    return acc + sum(t.weights .* (x .- r) .^ 2)
end

"""
    TerminalEllipsoidCost(E; w_outside = 1e6, w_center = 1e5)

Terminal cost pulling the endpoint into the ellipsoid `E` (typically the
certifier's terminal ellipsoid): with `d² = (x−c)ᵀP(x−c)`,
`w_outside·max(0, d²−1)² + w_center·d²`.
"""
struct TerminalEllipsoidCost{TC, TP} <: AbstractCostTerm
    c::TC
    P::TP
    w_outside::Float64
    w_center::Float64
end

function TerminalEllipsoidCost(E; w_outside = 1e6, w_center = 1e5)
    return TerminalEllipsoidCost(
        collect(Float64, LazySets.center(E)),
        Matrix{Float64}(UT.get_quadratic_form(E)),
        Float64(w_outside),
        Float64(w_center),
    )
end

function cost_final(t::TerminalEllipsoidCost, acc, xT)
    d = collect(xT) - t.c
    d2 = LA.dot(d, t.P * d)
    return acc + t.w_outside * max(0.0, d2 - 1.0)^2 + t.w_center * d2
end

"""
    InputEffortCost(w = 1.0)

Stage cost `w·Σₖ ‖uₖ‖²`.
"""
struct InputEffortCost <: AbstractCostTerm
    w::Float64
end
InputEffortCost() = InputEffortCost(1.0)

cost_step(t::InputEffortCost, acc, x, u, k) = acc + t.w * sum(abs2, u)

"""
    InputSmoothnessCost(w_du = 1.0, w_ddu = 0.0)

Stage cost on input increments: `w_du·‖Δuₖ‖² + w_ddu·‖Δ²uₖ‖²`. Smooth inputs are
what the certifier's linearization boxes can absorb (plan.md §3.2-M2).
"""
struct InputSmoothnessCost <: AbstractCostTerm
    w_du::Float64
    w_ddu::Float64
end
InputSmoothnessCost(; w_du = 1.0, w_ddu = 0.0) =
    InputSmoothnessCost(Float64(w_du), Float64(w_ddu))

cost_init(::InputSmoothnessCost) = (0.0, nothing, nothing)

function cost_step(t::InputSmoothnessCost, acc, x, u, k)
    total, u1, u2 = acc
    if u1 !== nothing
        du = u .- u1
        total += t.w_du * sum(abs2, du)
        if u2 !== nothing
            ddu = u .- 2 .* u1 .+ u2
            total += t.w_ddu * sum(abs2, ddu)
        end
    end
    return (total, u, u1)
end

cost_final(::InputSmoothnessCost, acc, xT) = acc[1]

"""
    DomainPenaltyCost(X; w = 1000.0, wrap = identity)

Stage penalty `w` per state outside the domain `X` (`wrap` maps periodic states
into the fundamental domain first). Every state of the rollout is charged —
including the terminal one — so leaving the domain early is never cheap.
"""
struct DomainPenaltyCost{TX, FW} <: AbstractCostTerm
    X::TX
    w::Float64
    wrap::FW
end
DomainPenaltyCost(X; w = 1000.0, wrap = identity) = DomainPenaltyCost(X, Float64(w), wrap)

cost_step(t::DomainPenaltyCost, acc, x, u, k) = acc + (t.wrap(x) ∈ t.X ? 0.0 : t.w)
cost_final(t::DomainPenaltyCost, acc, xT) = acc + (t.wrap(xT) ∈ t.X ? 0.0 : t.w)

"""
    ReachObjectiveCost(target_set; w_distance = 100.0, hit_bonus = 1000.0, wrap = identity)

Reach shaping: `w_distance · min_k dist(xₖ, target centers)` minus `hit_bonus / k`
for the first state inside the target (earlier hits earn more). Distances are to
the centers of the target's members — a cheap online surrogate.
"""
struct ReachObjectiveCost{TS, FW} <: AbstractCostTerm
    centers::Vector{Vector{Float64}}
    target_set::TS
    w_distance::Float64
    hit_bonus::Float64
    wrap::FW
end

function ReachObjectiveCost(
    target_set;
    w_distance = 100.0,
    hit_bonus = 1000.0,
    wrap = identity,
)
    members =
        target_set isa LazySets.UnionSetArray ? LazySets.array(target_set) : [target_set]
    centers = [Vector{Float64}(LazySets.center(m)) for m in members]
    return ReachObjectiveCost(
        centers,
        target_set,
        Float64(w_distance),
        Float64(hit_bonus),
        wrap,
    )
end

cost_init(::ReachObjectiveCost) = (Inf, 0)

function _reach_acc(t::ReachObjectiveCost, acc, x, k)
    best, hit = acc
    xw = t.wrap(x)
    d = minimum(LA.norm(collect(xw) - c) for c in t.centers)
    best = min(best, d)
    hit = (hit == 0 && xw ∈ t.target_set) ? k : hit
    return (best, hit)
end

cost_step(t::ReachObjectiveCost, acc, x, u, k) = _reach_acc(t, acc, x, k)

function cost_final(t::ReachObjectiveCost, acc, xT)
    best, hit = _reach_acc(t, acc, xT, typemax(Int))
    return t.w_distance * best - (hit == 0 ? 0.0 : t.hit_bonus / hit)
end

# Composite evaluation used by the rollout kernel.
_composite_init(c::CompositeCost) = map(cost_init, c.terms)
_composite_step(c::CompositeCost, accs, x, u, k) =
    map((t, a) -> cost_step(t, a, x, u, k), c.terms, accs)
_composite_final(c::CompositeCost, accs, xT) =
    sum(map((t, a) -> cost_final(t, a, xT), c.terms, accs))
