module DiscreteSystems

import Dionysos
const ST = Dionysos.System
const SY = Dionysos.Symbolic
const PR = Dionysos.Problem

import JuMP: MOI

import ..AbstractDionysosOptimizer

# The serializable q -> [u…] map lives in System next to DiscreteStaticController.
const DiscreteControlTable = ST.ControlTable
const add_control! = ST.add_control!
const set_control! = ST.set_control!

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

function _bitset_from_states(states, nstates::Int)
    bits = falses(nstates)
    for q in states
        bits[q] = true
    end
    return bits
end

function _set_from_bitset(bits::BitVector)
    out = Set{Int}()
    for q in eachindex(bits)
        bits[q] && push!(out, q)
    end
    return out
end

function _empty_control_table(nstates::Int)
    return DiscreteControlTable(nstates)
end

"""
    covers_initial_set(predicate, initial_set) -> Bool

Whether every abstract initial state satisfies `predicate` — the success criterion of a
control solver.

Unlike a bare `all`, an **empty** `initial_set` yields `false`. An empty set means the
specification's initial region is not represented in the abstraction at all (a degenerate
grid, a region outside the abstracted domain, …), so nothing was verified; reporting success
there is a vacuous truth that reaches the user as a spurious `MOI.OPTIMAL`.
"""
covers_initial_set(predicate, initial_set) =
    !isempty(initial_set) && all(predicate, initial_set)

# ------------------------------------------------------------
# Precompute existing (state,input) pairs once.
# base_pairstable[q,u] == true iff Post(q,u) is nonempty.
# ------------------------------------------------------------

function _compute_base_pairstable(autom)
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    base_pairstable = falses(nstates, nsymbols)

    for (_, source, symbol) in SY.enum_transitions(autom)
        base_pairstable[source, symbol] = true
    end

    return base_pairstable
end

function _compute_nsymbolslist(pairstable::BitMatrix)
    nstates, nsymbols = size(pairstable)
    nsymbolslist = zeros(Int, nstates)

    @inbounds for q in 1:nstates
        c = 0
        for u in 1:nsymbols
            c += pairstable[q, u]
        end
        nsymbolslist[q] = c
    end

    return nsymbolslist
end

function _counter_dense(autom)
    counter = zeros(Int, SY.get_n_state(autom), SY.get_n_input(autom))

    for target in SY.enum_states(autom)
        for (source, symbol) in SY.pre(autom, target)
            counter[source, symbol] += 1
        end
    end

    return counter
end

include("optimal_control_problem.jl")
include("safety_problem.jl")
include("reach_and_stay_problem.jl")
include("cosafe_ltl_problem.jl")

end
