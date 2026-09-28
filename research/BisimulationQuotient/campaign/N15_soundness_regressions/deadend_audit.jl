# Are the quotients total, and where do the deadend states sit?
#
# A bisimulation of a total system cannot contain a state with no outgoing transition, so a deadend
# is either a defect or an artefact of the `atol` inset. The distinction is settled by *where* the
# deadends sit and *how much volume* they carry: dust in the outermost slice is the inset leaving
# hairline residues, whereas a fat cell in an inner slice would be a defect.
#
# Measured result: three of the four cached quotients are total. The fourth has 3 deadends out of
# 9 027, carrying 3e-5 of a domain of volume 122.5 -- two ten-millionths of the state space, all in
# the unobserved region. That is inset dust, four orders of magnitude below the ~0.4 % of the domain
# the inset is already documented to discard at `atol = 1e-4`.
#
# It still has to be handled rather than left open: either prune cells below a volume threshold tied
# to `atol` so totality holds by construction, or state that the quotient is a bisimulation of the
# system restricted to the covered set, up to the measure the inset discards.
#
# This file also reports the slice count against the `max_slices` requested, which is the second
# open anomaly. It does not reproduce in any cached quotient.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

cell_volume(S::UT.SemiLinearSet) = sum(LazySets.volume, S.array)

"""
Deadend states of `quotient`, grouped by the slice they sit in, with the volume they carry.
"""
function deadend_report(quotient)
    dead = PCQ.deadend_states(quotient)
    n = length(quotient.states)
    println(
        "states = ",
        n,
        "   slices = ",
        PCQ.num_slices(quotient),
        "   deadends = ",
        length(dead),
        " (",
        round(100 * length(dead) / n; digits = 2),
        " %)",
    )
    if isempty(dead)
        println("  -> TOTAL: no deadend states.")
        return nothing
    end

    per_slice = Dict{Int, Int}()
    all_slices = Dict{Int, Int}()
    for state in values(quotient.states)
        all_slices[state.slice] = get(all_slices, state.slice, 0) + 1
    end
    for id in dead
        s = quotient.states[id].slice
        per_slice[s] = get(per_slice, s, 0) + 1
    end
    outermost = maximum(keys(all_slices))
    println("  deadends by slice (deadend / total in that slice):")
    for s in sort(collect(keys(per_slice)))
        tag = s == outermost ? "   <- outermost (truncation)" : "   <- inner"
        println("     slice ", s, " : ", per_slice[s], " / ", all_slices[s], tag)
    end

    dead_volume = sum(id -> cell_volume(quotient.states[id].set), dead)
    total_volume = sum(state -> cell_volume(state.set), values(quotient.states))
    println(
        "  deadend volume = ",
        dead_volume,
        " of ",
        round(total_volume; digits = 5),
        "  (",
        round(100 * dead_volume / total_volume; digits = 6),
        " %)",
    )
    return nothing
end

function report(name, path)
    if !isfile(path)
        println(
            "\n===== ",
            name,
            " =====\n  cache missing -- run the experiment that writes it",
        )
        return nothing
    end
    println("\n===== ", name, " =====")
    optimizer = import_optimizer_jld2(path)
    deadend_report(MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient")))
    return nothing
end

const ROOT = dirname(dirname(@__DIR__))

println("Totality audit of the cached quotients")
report("E2 PCLF", joinpath(ROOT, "experiments", "pclf_case.jld2"))
report("E2 induced CLF", joinpath(ROOT, "experiments", "clf_case.jld2"))
report("E8 PCLF", joinpath(ROOT, "examples", "gol_lazar_belta_pclf.jld2"))
report("E8 induced CLF", joinpath(ROOT, "examples", "gol_lazar_belta_clf.jld2"))
