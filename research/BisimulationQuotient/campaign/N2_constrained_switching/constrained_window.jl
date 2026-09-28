# Their matrices, made unstable: where does constrained switching become the only option?
#
# Scale Gol-Lazar-Belta's first mode, Σ_c = {c A₁, A₂}. For c above ≈ 1.17 the constant word 1^ω
# diverges, so the system is unstable under arbitrary switching and **no common Lyapunov function
# exists at any complexity** -- the predecessor's construction has no input at all. One eigenvalue
# establishes it, so unlike every search-based separation in this folder it cannot be defeated by
# giving the baseline a larger budget.
#
# Yet the language forbidding two consecutive uses of mode 1 remains stable well past that point, so
# a PCLF over the two-state automaton generating it still certifies. Same working set, same three
# observation regions, same specification: their problem, and their method cannot start.
#
# Measured window: c ∈ [1.17, 1.5].  Take c = 1.25, which has margin at both ends.
#
# What this file does NOT establish, and what must be checked before the experiment is built: the
# periodic bound says the constrained system is stabilizable at rate ≥ 0.912, not that a polyhedral
# template of a given order reaches it. Run the certificate search before committing.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

"""
Best periodic rate over words of length at most `Lmax`, optionally skipping words containing `11`.

`forbid_11` restricts the search to the constrained language: a cyclic word is admissible when no
two consecutive positions, wrapping around, both take mode 1.
"""
function best_periodic_rate(A, Lmax::Int; forbid_11::Bool = false)
    best, best_word = 0.0, Int[]
    for L in 1:Lmax, code in 0:(2 ^ L - 1)
        word = [((code >> (i - 1)) & 1) + 1 for i in 1:L]
        if forbid_11 && any(word[i] == 1 && word[mod1(i + 1, L)] == 1 for i in 1:L)
            continue
        end
        P = Matrix{Float64}(LinearAlgebra.I, 2, 2)
        for m in word
            P = A[m] * P
        end
        r = spectral_radius(P)^(1 / L)
        if r > best
            best, best_word = r, word
        end
    end
    return best, best_word
end

function verdict(arbitrary, constrained)
    if arbitrary < 1.0
        return "stable under arbitrary switching (uninteresting)"
    elseif constrained < 1.0
        return "*** unstable arbitrary, stable on no-11 ***"
    else
        return "unstable even on no-11"
    end
end

function sweep(cs; Lmax = 10)
    A1 = [-0.65 0.32; -0.42 -0.92]
    A2 = [0.65 0.32; -0.42 -0.92]
    println(
        "Σ_c = {c A₁, A₂} with Gol-Lazar-Belta's matrices; periodic bounds over L ≤ ",
        Lmax,
    )
    println(
        rpad("c", 7),
        rpad("ρ(cA₁)", 12),
        rpad("JSR lb arbitrary", 20),
        rpad("JSR lb no-11", 18),
        "verdict",
    )
    for c in cs
        A = [c * A1, A2]
        arbitrary, _ = best_periodic_rate(A, Lmax)
        constrained, _ = best_periodic_rate(A, Lmax; forbid_11 = true)
        println(
            rpad(c, 7),
            rpad(round(spectral_radius(c * A1); digits = 4), 12),
            rpad(round(arbitrary; digits = 4), 20),
            rpad(round(constrained; digits = 4), 18),
            verdict(arbitrary, constrained),
        )
    end
    return nothing
end

sweep([1.0, 1.1, 1.17, 1.25, 1.3, 1.4, 1.5, 1.6])
