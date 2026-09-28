# Why refining a uniform grid cannot reduce its nondeterminism.
#
# For a linear map and a tight box over-approximation, a cell of radius h/2 maps to a box of radius
# |A|h/2, which needs about prod_i(sum_j |A_ij| + 1) cells to cover -- the image and the cell scale
# together, so the branching factor is INDEPENDENT of h.
#
# Measured on Gol-Lazar-Belta, mean successors per (state, mode):
#
#   h = 1.0    377 cells      4.59        h = 0.25    6 591 cells    4.61
#   h = 0.5  1 603 cells      4.62        h = 0.125  26 751 cells    4.61
#
# Flat at 4.61 across a 70x increase in cell count, matching the closed form 1.97 * 2.34 = 4.61 to
# three digits, with no cell losing all its transitions. Refining buys precision in WHERE states are
# and none at all in HOW NONDETERMINISTIC the abstraction is: the adversary keeps ~4.6 choices per
# step forever.
#
# Contrast the bisimulation quotients, measured at 1 successor per (state, mode) whenever the graph
# is deterministic (0 of 21 220 pairs branch). A quotient branches only when its GRAPH does, and that
# branching is real structure rather than error.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))
import MathematicalSystems
const SY = DI.Symbolic

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
A = UT.mode_matrices(f)
mode_of(u) = clamp(round(Int, u[1]), 1, length(A))
step_map(x, u) = A[mode_of(u)] * x
function img(rect, u)
    c = LazySets.center(rect)
    r = LazySets.radius_hyperrectangle(rect)
    M = A[mode_of(u)]
    return LazySets.Hyperrectangle(M * c, abs.(M) * r)
end

CACHE = joinpath(dirname(dirname(@__DIR__)), "examples", "gol_lazar_belta_pclf.jld2")
opt0 = import_optimizer_jld2(CACHE)
sets = [PCLF.get_sublevel_set(p, maximum(opt0.Γ)) for (_, p) in opt0.pclf.pieces]
domain = sets[argmax([LazySets.volume(s) for s in sets])]

for (i, M) in enumerate(A)
    rs = sum(abs.(M); dims = 2)
    println(
        "mode ",
        i,
        ": |A| row sums = ",
        round.(vec(rs); digits = 3),
        "  -> predicted successors ~ ",
        round(prod(vec(rs) .+ 1); digits = 2),
    )
end
println()

function branch(h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        domain,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    o = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        o,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(system, domain),
    )
    MOI.set(
        o,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0, 0.0), SVector(h, h)),
    )
    MOI.set(
        o,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(1.0), SVector(1.0)),
    )
    MOI.set(
        o,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(o, MOI.RawOptimizerAttribute("overapproximation_map"), img)
    MOI.set(o, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(o)
    abs_sys = MOI.get(o, MOI.RawOptimizerAttribute("abstract_system"))
    autom = SY.get_automaton(abs_sys)
    n = SY.get_n_state(abs_sys)
    nt = SY.ntransitions(autom)
    # successors per (state, input): count transitions grouped by (source, symbol)
    per = Dict{Tuple{Int, Int}, Int}()
    for (tgt, src, sym) in SY.enum_transitions(autom)
        per[(src, sym)] = get(per, (src, sym), 0) + 1
    end
    vals = collect(values(per))
    dead = n - length(unique(first.(keys(per))))
    return n,
    nt,
    isempty(vals) ? 0.0 : sum(vals)/length(vals),
    isempty(vals) ? 0 : maximum(vals),
    length(per),
    dead
end

println(
    rpad("h", 9),
    rpad("cells", 9),
    rpad("transitions", 13),
    rpad("mean succ", 11),
    rpad("max succ", 10),
    rpad("(s,u) pairs", 13),
    "cells w/o any transition",
)
for h in [1.0, 0.5, 0.25, 0.125]
    n, nt, mean_s, max_s, pairs, dead = branch(h)
    println(
        rpad(h, 9),
        rpad(n, 9),
        rpad(nt, 13),
        rpad(round(mean_s; digits = 2), 11),
        rpad(max_s, 10),
        rpad(pairs, 13),
        dead,
    )
end
