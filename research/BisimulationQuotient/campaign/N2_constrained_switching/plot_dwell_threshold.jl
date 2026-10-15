# Figure A for N2: how much conservatism buys a certificate.
#
# The narrative in one plot. The plant is unstable under arbitrary switching, so the predecessor's
# construction has no input at all -- there is no common Lyapunov function at any complexity, and the
# unconstrained rate sits above 1 however rich the template. Restrict the switching signal with a
# minimum dwell time and the rate falls through 1 at tau = 3: from that point on a certificate
# exists, a finite bisimulation can be built, and co-safe LTL synthesis is available.
#
# The second curve is the price. The dwell automaton has M*tau nodes, so the lifted state space
# carries that many copies of the domain. The figure therefore shows both halves of the trade at
# once: the certificate improves with tau, and the abstraction grows with tau.
#
# tau = 1 is the control: the constraint is vacuous there, and the rate must reproduce the
# unconstrained one. It does, to six digits.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

include(joinpath(@__DIR__, "dwell_time_example.jl"))

gr()

const TAUS = 1:6
const ORDER = 2

A1, A2 = shear_modes(; diag = 0.3, shear = 2.0)
f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([A1, A2]),
    HybridSystems.ControlledSwitching(),
)

unconstrained_rate = rate(f, PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 1, 2)]), ORDER)
rates = [rate(f, dwell_graph(2, τ), ORDER) for τ in TAUS]
nodes = [2 * τ for τ in TAUS]

println("unconstrained: ", round(unconstrained_rate; digits = 5))
for (τ, r) in zip(TAUS, rates)
    println("  tau = ", τ, "  nodes = ", 2τ, "  rate = ", round(r; digits = 5))
end

fig = plot(;
    xlabel = "minimum dwell time τ",
    ylabel = "certified contraction rate",
    legend = :topright,
    size = (620, 420),
    title = "Conservatism buys a certificate",
)

# The unstable region: any rate at or above 1 certifies nothing.
hspan!(
    fig,
    [1.0, maximum([unconstrained_rate; rates]) * 1.05];
    color = LOSING_COLOR,
    alpha = 0.10,
    label = "no certificate",
)
hline!(fig, [1.0]; color = :black, linestyle = :dash, linewidth = 1.5, label = "")

hline!(
    fig,
    [unconstrained_rate];
    color = LOSING_COLOR,
    linewidth = 2,
    linestyle = :dot,
    label = "unconstrained (no common Lyapunov function exists)",
)
plot!(
    fig,
    collect(TAUS),
    rates;
    marker = :circle,
    linewidth = 2,
    color = WINNING_COLOR,
    label = "dwell-time constrained",
)

# Annotate the threshold: the first tau whose rate certifies.
τ_star = findfirst(<(1.0), rates)
if τ_star !== nothing
    scatter!(
        fig,
        [TAUS[τ_star]],
        [rates[τ_star]];
        marker = (:star5, 12),
        color = WINNING_COLOR,
        label = "τ* = $(TAUS[τ_star]) ($(nodes[τ_star]) nodes)",
    )
end

savefig(fig, joinpath(@__DIR__, "fig_dwell_threshold.png"))
println("\nwrote fig_dwell_threshold.png")
