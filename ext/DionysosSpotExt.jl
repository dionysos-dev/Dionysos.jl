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

# The good-prefix states of the monitor: those from which *every* infinite extension is accepted,
# which for a co-safe formula is exactly "the prefix so far already satisfies φ".
#
# This replaces an absorbing-states heuristic that read a missing edge as a self-loop. A missing
# edge is a step into the rejecting dead state, so that heuristic could declare a state accepting
# although some extension from it violates φ. Under synthesis the resulting controller fails
# visibly; under verification (a folded run) the mistake silently verifies violating states —
# which is why this must be computed, not guessed.
#
# The computation is a pair of greatest fixed points over the explicit valuation alphabet:
#
#   1. `C` — the largest set of reachable states that is *total* (every valuation has an edge)
#      and closed (every edge stays in `C`). A run entering `C` never dies.
#   2. For each Rabin pair `(Fin, Inf)`: the largest closed subset of `C ∩ (Inf ∖ Fin)`. Every
#      infinite run inside such a subset visits `Inf` forever and `Fin` never, so it is accepted
#      whatever the word — the certificate that every extension is good.
#
# The union over pairs under-approximates the true good-prefix set in general and is exact on the
# monitors spot produces for co-safe formulas, whose good prefixes end in a terminal accepting
# component. Under-approximation is the sound direction for both quantifiers: a smaller accepting
# set can only lose controllers or verified states, never invent them.
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

    # Greatest fixed point: keep the states whose every valuation stays inside the kept set.
    function largest_closed_subset(candidates::Set{Int})
        S = copy(candidates)
        changed = true
        while changed
            changed = false
            for q in collect(S)
                for v in vals
                    q2 = Spot.nextstate(dra, q, v)
                    if q2 === nothing || !(q2 in S)
                        delete!(S, q)
                        changed = true
                        break
                    end
                end
            end
        end
        return S
    end

    C = largest_closed_subset(reachable)

    good = Set{Int}()
    for (fin, inf) in Spot.get_rabin_acceptance(dra)
        candidates = Set(q for q in C if (q in inf) && !(q in fin))
        union!(good, largest_closed_subset(candidates))
    end

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
    goodQ, _ = _good_prefix_states_dra(S.dra, vals; q0 = S.qa0)
    isempty(goodQ) && error(
        "No state of the monitor certifies a good prefix over the emitted labels: no closed, " *
        "total, always-accepting component is reachable. Either the formula has no good prefix " *
        "over this alphabet, or the monitor's acceptance is beyond the terminal certificate " *
        "computed here — provide a FunctionMonitor with explicit accepting states rather than " *
        "guessing.",
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
