# Experiment 2 for N10, on Gol-Lazar-Belta's own geometry and their own co-safe LTL formula.
#
# The dwell-time example of `quotient_on_dwell.jl` makes the same point but renders badly: its shear
# makes the certificate's sublevel sets wildly eccentric -- they reach about ±130 against a ±4 working
# set -- so it is kept as a structural diagnostic with no figure. Their matrices are well proportioned,
# which is why the E8 figures read well, so this folder's figure is built on those instead.
#
# Everything is theirs except one scalar and the constraint:
#
#   system   {c*A1, A2} with c = 1.20.  rho(c*A1) = 1.0270 > 1, so the constant word 1^omega
#            diverges and NO common Lyapunov function exists at any complexity -- their construction
#            has no input at all.  One eigenvalue establishes it; no search can rescue it.
#   language no two consecutive uses of mode 1, a two-state automaton.  Certifies at 0.901.
#   working set, three observation regions, and specification: exactly theirs.
#
# c is tuned, and the trade is worth recording: a smaller c gives a better constrained rate, hence
# fewer slices and a smaller quotient, at the cost of a thinner instability margin.
#
#   c      rho(c*A1)   margin   no-11 rate (order 2)
#   1.18   1.0098      1.0 %    0.894289
#   1.20   1.0270      2.7 %    0.901091      <- used here
#   1.22   1.0441      4.4 %    0.907839
#   1.25   1.0698      7.0 %    0.917851      <- first attempt; quotient passed 25 000 cells
#
# Note the margin is a PROOF, not a measurement: rho(1.2*A1) = 1.027 > 1 means the word 1^omega
# diverges, however thin the margin looks.
#
#       phi = (!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D)
#
# So the figure is: their problem, their formula, and their method cannot start -- while restricting
# the switching signal yields a certificate, a finite bisimulation, and a certified region.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Spot
using LinearAlgebra

const C = 1.20
const ORDER = 2

(; problem, X, R1, R2, R3) = gol_lazar_belta_problem()
A1, A2 = UT.mode_matrices(problem.system)

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

# The scaled system: mode 1 alone now diverges.
scaled = ST.with_switching(
    HybridSystems.discreteswitchedsystem([C .* Matrix(A1), Matrix(A2)]),
    HybridSystems.ControlledSwitching(),
)
constrained_problem = PR.BisimulationQuotientProblem(scaled, X, [R1, R2, R3])

println(
    "rho(c*A1) = ",
    round(spectral_radius(C .* A1); digits = 4),
    " > 1  ->  no common Lyapunov function exists at any complexity",
)

# The language: no two consecutive 1s.  q0 --1--> q1, q0 --2--> q0, q1 --2--> q0.
no_11 = PCLF.edgeList_to_LabDigraph([(0, 1, 1), (0, 0, 2), (1, 0, 2)])

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "verbose" => false,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

pclf = PCLF.compute_polyhedral_pieces_pclf(
    scaled,
    no_11,
    OPTIMIZER,
    Dict(v => PCLF.conic_partitions_2d(ORDER) for v in no_11.verts);
    MLF = true,
)
println(
    "certificate on the constrained language: rate = ",
    round(pclf.JSRapprox; digits = 6),
)
isfinite(pclf.JSRapprox) || error("no certificate; retune c or the template order")

# Cached, so redrawing the figure does not pay for the build again. Delete the `.jld2` after
# changing `C`, `ORDER`, the graph or the tolerances, or a stale quotient is loaded against a fresh
# certificate.
const CACHE = joinpath(@__DIR__, "glb_constrained_c$(C).jld2")

quotient, D = if isfile(CACHE)
    println("reusing the cached quotient from ", basename(CACHE))
    opt = import_optimizer_jld2(CACHE)
    MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
    MOI.get(opt, MOI.RawOptimizerAttribute("D"))
else
    r = build_quotient(
        constrained_problem,
        pclf;
        atol = 1e-4,
        level_tol = 1e-2,
        max_slices = 40,
        print_level = 1,
    )
    export_optimizer_jld2(r.optimizer, CACHE)
    r.quotient, r.D
end
println(
    "\nquotient: ",
    length(quotient.states),
    " cells on ",
    length(no_11.verts),
    " nodes",
)

# Their formula, unchanged.
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"

# Their initial point `a`, exactly as `examples/gol_lazar_belta_pclf.jl` uses it.
const X0 = SVector(-4.0, -7.0)

"""
The certified set, without needing a valid initial point.

`synthesize_cosafe_ltl` also builds a controller and simulates it, which requires an initial state
that is both covered by the quotient and winning. Neither is guaranteed here: `max_slices` truncates
the covered region, so their `a = (-4, -7)` can fall outside it, and a point inside the terminal set
is covered but cannot discharge the formula's `F(R1)` obligation. The figure needs only the certified
*set*, so this runs the synthesis and stops before the controller is built.
"""
function certified_set(quotient, spec, regions, ap_to_obs)
    # A POINT, not the working set. The co-safe product's initial MEMORY state is fixed by the
    # observation at the initial set, and `project_initial_memory_controllable_set` reports the
    # certified set at that memory. `X` spans all three regions at once, so it starts the monitor
    # from the wrong state -- silently, as a merely small answer. This script passed `X` until the
    # figure was redrawn: 53 of 2 138 cells. `glb_constraint_moves_both_ways.jl` hit the same thing
    # on the unconstrained control, where it showed as 442 against E8's known 8 794.
    initial_set = LazySets.Hyperrectangle(; low = [X0[1], X0[2]], high = [X0[1], X0[2]])
    problem = PR.CoSafeLTLProblem(
        scaled,
        initial_set,
        spec,
        regions,
        Dict{Symbol, Any}(ap => MP.INNER for ap in keys(regions)),
    )
    optimizer = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"), quotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("ap_to_obs"), ap_to_obs)
    # `early_stop = false` is what makes this the certified set of the WHOLE domain rather than of
    # the initial set alone -- the optimizer's default is `true`, which seeds the product only from
    # states meeting the initial set. It matters more here than anywhere: the covered region reaches
    # ±17 while `X` is ±10, so with `true` every cell outside the working set is unreachable from the
    # seed and can never be certified, however good the certificate is.
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    return MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
end

win = collect(
    certified_set(
        quotient,
        Dionysos.spot_stepper(φ),
        Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
        Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    ),
)
lose = collect(setdiff(Set(keys(quotient.states)), Set(win)))
println(
    "certified for their formula: ",
    length(win),
    " of ",
    length(quotient.states),
    " cells",
)

gr()

# No fixed window. The certificate's sublevel sets reach well past the ±10 working set -- the covered
# region is a star, not a box -- so the hard-coded ±10 this script used to carry showed only the
# middle of the result. The quotient panel auto-scales and the other panel is matched to it, so the
# two stay comparable. (`quotient_on_dwell.jl` does clip on purpose: its shear pushes the sublevel
# sets out to about ±130 against a ±4 working set, where auto-scaling shows nothing legible.)
ours = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "constrained: $(length(win))/$(length(quotient.states)) certified",
    titlefontsize = 10,
)
# Contours ON, which looks backwards and is not. The white haze `linewidth = 0` leaves over a panel
# this dense is not strokes -- the recipe already forces `linealpha = 0` when `show_contours` is
# false -- it is GR's anti-aliased seam between each pair of adjacent filled polygons. Turning
# contours on strokes every cell in `linecolor := c`, its own fill colour, covering the seams.
#
# The width has to beat the seam at THIS panel's scale, and the scale differs per figure: 0.7 clears
# `fig_constraint_both_semantics.png`, whose cells span a ±13 window, and leaves a visible haze here,
# where the same machinery draws a ±17 window into a panel of the same size.
_plot_winning_losing!(
    ours,
    quotient,
    win,
    lose,
    nothing;
    show_contours = true,
    linewidth = 1.5,
)
plot!(
    ours,
    constrained_problem;
    plot_region = false,
    observation_region_alpha = 0.0,
    observation_colors = OBSERVATION_COLORS,
    observation_linewidth = 2.0,
)
LIMS_X = Plots.xlims(ours)
LIMS_Y = Plots.ylims(ours)
println(
    "panel window: x ∈ ",
    round.(LIMS_X; digits = 1),
    ", y ∈ ",
    round.(LIMS_Y; digits = 1),
)

# Their panel is empty by construction, and the annotation says exactly why.
theirs = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "Gol–Lazar–Belta: no construction",
    titlefontsize = 10,
    xlims = LIMS_X,
    ylims = LIMS_Y,
)
plot!(theirs, X; color = :grey90, linecolor = :black, fillalpha = 1.0)
plot!(
    theirs,
    constrained_problem;
    plot_region = false,
    observation_region_alpha = 0.0,
    observation_colors = OBSERVATION_COLORS,
    observation_linewidth = 2.0,
)
# Below the working set, not on top of it. Sharing the quotient panel's window shrinks X to the
# middle of this one, so an annotation centred at the origin lands squarely on the three region
# outlines and is unreadable there. The lower third of the panel is empty, which is the point.
annotate!(
    theirs,
    0.0,
    LIMS_Y[1] + 0.09 * (LIMS_Y[2] - LIMS_Y[1]),
    text(
        "ρ(1.20·A₁) = 1.027 > 1\nmode 1 alone diverges\nno common Lyapunov function\nat any complexity",
        9,
        :center,
    ),
)

fig = plot(theirs, ours; layout = (1, 2), size = (960, 480))
savefig(fig, joinpath(@__DIR__, "fig_glb_constrained.png"))
println("\nwrote fig_glb_constrained.png")
