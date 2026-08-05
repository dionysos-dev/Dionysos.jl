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

function _cosafe_done_states_dra(dra; aps::Vector{Symbol}, q0::Int = 1)
    vals = _all_valuations(aps)

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

    done = Set{Int}()
    for q in reachable
        ok = true
        for v in vals
            q2 = Spot.nextstate(dra, q, v)
            if !(q2 === nothing || q2 == q)
                ok = false
                break
            end
        end
        ok && push!(done, q)
    end

    return done, reachable
end

OPDS.accepting_states(S::SpotDRAstepper) = begin
    doneQ, _ = _cosafe_done_states_dra(S.dra; aps = S.aps, q0 = S.qa0)
    isempty(doneQ) && error(
        "Could not find any 'done' states with the co-safe heuristic. " *
        "Try providing a FunctionMonitor with explicit accepting states.",
    )
    doneQ
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
