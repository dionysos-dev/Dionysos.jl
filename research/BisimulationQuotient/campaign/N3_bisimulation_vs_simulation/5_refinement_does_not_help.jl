# Does refining the grid eventually let it satisfy the co-safe formula?  No -- it gets worse.
#
#   h = 0.25     6 591 cells    48 certified   area 3.00   0.71 % of the domain
#   h = 0.125   26 751 cells    83 certified   area 1.30   0.31 %
#   h = 0.0625 107 791 cells   100 certified   area 0.39   0.09 %
#
# At ten times the bisimulation's 10 611 cells the grid manages 0.09 %, against the bisimulation's
# 83 %. The certified AREA decays roughly in proportion to h even as the certified CELL count rises,
# which is the signature of a boundary layer: the fixed point propagates about one step out from the
# accepting set and no further.
#
# The script also prints why resolution is not the binding constraint: the terminal set D has area
# 141.3 of 424 -- a third of the domain, equivalent radius 6.7 -- while the residual quantization
# uncertainty h/(2(1-rho)) = 3.77h is only 0.24 at the finest step. The target is enormous and the
# grid still cannot reach it.  See `4_branching_is_scale_invariant.jl` for the reason.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))
using Spot
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
D = MOI.get(opt0, MOI.RawOptimizerAttribute("D"))
sets = [PCLF.get_sublevel_set(p, maximum(opt0.Γ)) for (_, p) in opt0.pclf.pieces]
domain = sets[argmax([LazySets.volume(s) for s in sets])]
regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"

Darea = sum(LazySets.volume, D.array)
println(
    "terminal set D: ",
    length(D.array),
    " parts, area ",
    round(Darea; digits = 3),
    "  (equivalent disc radius ",
    round(sqrt(Darea/pi); digits = 3),
    ")",
)
println(
    "rho = ",
    opt0.pclf.JSRapprox,
    " -> residual uncertainty ~ h/(2(1-rho)) = ",
    round(1/(2*(1-opt0.pclf.JSRapprox)); digits = 2),
    " * h\n",
)

function run(h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        domain,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    concrete = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = [-4.0, -7.0], high = [-4.0, -7.0]),
        Dionysos.spot_stepper(φ),
        regions,
        Dict{Symbol, Any}(:D=>MP.INNER, :R1=>MP.INNER, :R2=>MP.OUTER, :R3=>MP.OUTER),
    )
    o = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(o, MOI.RawOptimizerAttribute("concrete_problem"), concrete)
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
    t = @elapsed MOI.optimize!(o)
    abs_sys = MOI.get(o, MOI.RawOptimizerAttribute("abstract_system"))
    n = SY.get_n_state(abs_sys)
    w = length(o.control_solver.abstract_optimizer.controllable_set)
    return n, w, w*h^2, t
end

println(
    rpad("h", 9),
    rpad("cells", 9),
    rpad("certified", 11),
    rpad("area", 9),
    rpad("% of 424", 10),
    "time (s)",
)
for h in [0.25, 0.125, 0.0625]
    n, w, a, t = run(h)
    println(
        rpad(h, 9),
        rpad(n, 9),
        rpad(w, 11),
        rpad(round(a; digits = 2), 9),
        rpad(round(100*a/424.12; digits = 2), 10),
        round(t; digits = 1),
    )
end
