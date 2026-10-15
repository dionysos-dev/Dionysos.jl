# Does the grid fail on Gol-Lazar-Belta because of the DYNAMICS or because of the SPECIFICATION?
#
# Controlled comparison: same matrices, same domain, same grid, same 4.61-way branching -- only the
# specification changes. Plain reachability to D and to R1, against the co-safe formula of
# `3_cosafe_glb.jl`.
#
# Result: it is the specification, decisively.
#
#   reach D    h=0.5   94.5 %    h=0.25  97.1 %    h=0.125  98.6 %   of the domain
#   reach R1   h=0.5   23.6 %    h=0.25  36.7 %    h=0.125  43.7 %   monotone, converging from below
#   their phi  h=0.25   0.71 %   h=0.125  0.31 %                     and falling
#
# So 4.61-way spurious branching is perfectly survivable for reachability. What it destroys is the
# ORDER of events, which is all the co-safe formula is about.
#
# Also prints the branching constant prod_i(sum_j |A_ij| + 1) for both benchmark systems: 2.97/2.89
# for two_mode, 4.61 for Gol-Lazar-Belta.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))
import MathematicalSystems
const SY = DI.Symbolic
const OPDS = OP.DiscreteSystems

for (nm, pr) in
    [("two_mode", two_mode_problem()), ("gol_lazar_belta", gol_lazar_belta_problem())]
    for (i, M) in enumerate(UT.mode_matrices(pr.f))
        rs = vec(sum(abs.(M); dims = 2))
        println(
            rpad(nm, 18),
            " mode ",
            i,
            ": row sums ",
            round.(rs; digits = 3),
            " -> branching ~ ",
            round(prod(rs .+ 1); digits = 2),
        )
    end
end
println()

# GLB matrices, PLAIN REACHABILITY to D: same dynamics, simpler specification.
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
D = MOI.get(opt0, MOI.RawOptimizerAttribute("D"))
sets = [PCLF.get_sublevel_set(p, maximum(opt0.Γ)) for (_, p) in opt0.pclf.pieces]
domain = sets[argmax([LazySets.volume(s) for s in sets])]

function reach(target, h)
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
    mp = SY.get_state_mapping(abs_sys)
    tg = collect(MP.get_states_from_set(mp, target, MP.INNER))
    _, ctrl, _, _ = OPDS.compute_worst_case_cost_controller(autom, tg)
    return SY.get_n_state(abs_sys), length(tg), length(ctrl), length(ctrl)*h^2
end

println("GLB matrices, PLAIN reachability (no monitor). Domain area 424.12.")
println(
    rpad("target", 9),
    rpad("h", 9),
    rpad("cells", 9),
    rpad("target", 9),
    rpad("certified", 11),
    rpad("area", 9),
    "% domain",
)
for (tn, tgt) in [("D", D), ("R1", R1)], h in [0.5, 0.25, 0.125]
    n, nt, c, a = reach(tgt, h)
    println(
        rpad(tn, 9),
        rpad(h, 9),
        rpad(n, 9),
        rpad(nt, 9),
        rpad(c, 11),
        rpad(round(a; digits = 2), 9),
        round(100*a/424.12; digits = 2),
    )
end
