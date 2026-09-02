module DionysosSpotExt

import Dionysos
import Spot
import Dionysos: spot_stepper

const DI = Dionysos
const OPDS = DI.Optim.DiscreteSystems

# ------------------------------------------------------------
# Spot helpers
# ------------------------------------------------------------

collect_aps(φ::Spot.SpotFormula) = [Symbol(ap) for ap in Spot.atomic_prop_collect(φ)]

struct SpotDRAstepper <: OPDS.AbstractSpecStepper
    φ::Spot.SpotFormula
    dra::Any
    qa0::Int
    qa_dead::Int
    aps::Vector{Symbol}
end

OPDS.init_state(S::SpotDRAstepper) = S.qa0

@inline function _nextstate_int(dra, qa::Int, ap::Tuple, qa_dead::Int)
    qa2 = Spot.nextstate(dra, qa, ap)
    return qa2 === nothing ? qa_dead : qa2
end

OPDS.step(S::SpotDRAstepper, qa::Int, ap::Tuple{Vararg{Symbol}}) =
    (qa == S.qa_dead) ? S.qa_dead : _nextstate_int(S.dra, qa, ap, S.qa_dead)

function _all_valuations(aps::Vector{Symbol})
    n = length(aps)
    vals = Vector{Tuple{Vararg{Symbol}}}(undef, 1 << n)
    for mask in 0:((1 << n) - 1)
        trues = Symbol[]
        for (i, a) in enumerate(aps)
            if (mask >> (i - 1)) & 1 == 1
                push!(trues, a)
            end
        end
        vals[mask + 1] = Tuple(trues)
    end
    return vals
end

# The accepting states of the monitor under the FINITE-TRACE reading the solvers implement: a
# closed loop stops the moment it reaches acceptance, so a state accepts iff the run stopped
# there — padded with the empty valuation forever, since nothing further is emitted — is an
# accepted word. Obligations (`F`, `U`) must already be discharged, prohibitions (`G(!p)`) hold
# vacuously over the padding.
#
# This replaces an absorbing-states heuristic that read a missing edge as a self-loop. A missing
# edge is a step into the rejecting dead state, so that heuristic could declare a state accepting
# although the prefix does not satisfy the formula. Under synthesis the resulting controller
# fails visibly; under verification (a folded run) the mistake silently verifies violating
# states — which is why this must be computed, not guessed.
#
# On genuinely co-safe formulas this coincides with the good-prefix set (the padding lands in
# the monitor's terminal accepting component); its added reach is the reach-avoid shape
# `F(goal) & G(!hazard)`, whose "goal seen, hazard never" state accepts under the stop reading
# although no finite prefix certifies the infinite `G`. The caller guards the degenerate end of
# that spectrum: a formula whose EMPTY continuation is already accepted from the initial state
# is pure safety, has no finite-trace content, and is refused (`_checked_good_prefix`).
#
# The empty-valuation lasso from `q`: accepted iff its cycle meets some Rabin pair.
function _pad_lasso_accepted(dra, q0pad::Int)
    seen = Int[]
    q = q0pad
    while !(q in seen)
        push!(seen, q)
        q2 = Spot.nextstate(dra, q, ())
        q2 === nothing && return false
        q = q2
    end
    cycle = seen[findfirst(==(q), seen):end]
    for (fin, inf) in Spot.get_rabin_acceptance(dra)
        any(qc -> qc in fin, cycle) && continue
        any(qc -> qc in inf, cycle) && return true
    end
    return false
end

function _good_prefix_states_dra(dra, vals; q0::Int = 1)
    reachable = Set{Int}([q0])
    queue = [q0]
    while !isempty(queue)
        q = popfirst!(queue)
        for v in vals
            q2 = Spot.nextstate(dra, q, v)
            q2 === nothing && continue
            if !(q2 in reachable)
                push!(reachable, q2)
                push!(queue, q2)
            end
        end
    end

    good = Set{Int}(q for q in reachable if _pad_lasso_accepted(dra, q))
    return good, reachable
end

# Without alphabet information the free powerset of the formula's atomic propositions is the only
# safe assumption; the product passes the labels the system can actually emit, which is both the
# correct notion of "every extension" and what keeps mutually exclusive observations from being
# judged against conjunctions that can never occur.
OPDS.accepting_states(S::SpotDRAstepper) = _checked_good_prefix(S, _all_valuations(S.aps))

OPDS.accepting_states(S::SpotDRAstepper, used_labels) =
    _checked_good_prefix(S, collect(Tuple{Vararg{Symbol}}, used_labels))

function _checked_good_prefix(S::SpotDRAstepper, vals)
    _pad_lasso_accepted(S.dra, S.qa0) && error(
        "The empty continuation already satisfies the formula from its initial state: this is " *
        "a safety property with no finite-trace content, and declaring everything accepting " *
        "would be unsound under verification. State it as a `SafetyProblem` (or add a " *
        "reachability obligation) instead.",
    )
    goodQ, _ = _good_prefix_states_dra(S.dra, vals; q0 = S.qa0)
    isempty(goodQ) && error(
        "No state of the monitor accepts under the finite-trace reading over the emitted " *
        "labels: no reachable state has its obligations discharged. Either the formula's " *
        "obligations cannot be met over this alphabet, or its acceptance is beyond the " *
        "stop-and-pad reduction computed here — provide a FunctionMonitor with explicit " *
        "accepting states rather than guessing.",
    )
    return goodQ
end

function spot_stepper(
    φ::Spot.SpotFormula;
    qa0::Int = 1,
    qa_dead::Int = 0,
    aps::Union{Nothing, Vector{Symbol}} = nothing,
)
    aps_use = aps === nothing ? collect_aps(φ) : aps
    dra = Spot.DeterministicRabinAutomata(φ)
    return SpotDRAstepper(φ, dra, qa0, qa_dead, aps_use)
end

# Lets the JuMP front-end take a formula directly: `@specification(model, ltl"F(goal)")`.
DI.Wrapper.to_stepper(φ::Spot.SpotFormula) = spot_stepper(φ)

end
