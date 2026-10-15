# The same comparison on Gol-Lazar-Belta's own example, with their co-safe LTL formula.
#
# `grid_vs_quotient.jl` compares the two abstractions on plain reachability, which is the one
# specification where the grid is least disadvantaged: no monitor, no product, and the fixed point
# is the same on both sides. This file removes that objection by asking their actual question,
#
#     φ = (!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D)
#
# on their system, their working set and their three observation regions.
#
# Note what D is. It is the *terminal set of the certificate*, not a region of the problem, so the
# grid has no natural counterpart for it -- the formula is only meaningful once a certificate has
# supplied one. The PCLF's D is therefore handed to the grid as a labelled region, which is the
# fairest reading available: both sides answer the same formula over the same atomic propositions.
# It is also worth saying plainly in the paper, because it is an argument in the method's favour
# that the specification itself is phrased in terms the certificate provides.
#
# Both abstractions discretize the same set: the certificate's outermost sublevel set.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Spot
using LinearAlgebra
import MathematicalSystems

const SY = DI.Symbolic
const OPDS = OP.DiscreteSystems

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const A = UT.mode_matrices(f)

mode_of(u) = clamp(round(Int, u[1]), 1, length(A))
step_map(x, u) = A[mode_of(u)] * x

function exact_image_box(rect, u)
    c = LazySets.center(rect)
    r = LazySets.radius_hyperrectangle(rect)
    M = A[mode_of(u)]
    return LazySets.Hyperrectangle(M * c, abs.(M) * r)
end

# ---------------------------------------------------------
# The quotient: E8's certificate, from cache when it exists
# ---------------------------------------------------------

const CACHE = joinpath(dirname(dirname(@__DIR__)), "examples", "gol_lazar_belta_pclf.jld2")

function load_or_build()
    if isfile(CACHE)
        println("reusing the cached E8 quotient")
        opt = import_optimizer_jld2(CACHE)
        return opt,
        MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
        MOI.get(opt, MOI.RawOptimizerAttribute("D"))
    end
    println("building the E8 quotient (a few minutes) ...")
    pclf = PCLF.compute_polyhedral_pieces_pclf(
        f,
        PCLF.generate_DeBruijn_edges(2, 1),
        JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
        Dict((1,) => PCLF.conic_partitions_2d(2), (2,) => PCLF.conic_partitions_2d(2));
        MLF = true,
    )
    r = build_quotient(problem, pclf; atol = 1e-4, level_tol = 1e-2, max_slices = 50)
    return r.optimizer, r.quotient, r.D
end

optimizer, quotient, D = load_or_build()

"""
The region the quotient covers: the largest of the per-node outermost sublevel sets.

The E8 pieces are nested, so their union is the outer one; that is asserted rather than assumed.
"""
function certificate_domain(pclf, Γ)
    sets = [PCLF.get_sublevel_set(p, maximum(Γ)) for (_, p) in pclf.pieces]
    volumes = [LazySets.volume(s) for s in sets]
    outer = sets[argmax(volumes)]
    for s in sets
        LazySets.issubset(s, outer) ||
            @warn "pieces are not nested; the union is larger than the biggest sublevel set"
    end
    return outer
end

domain = certificate_domain(optimizer.pclf, optimizer.Γ)

# ---------------------------------------------------------
# The specification: theirs
# ---------------------------------------------------------

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"

# Polarity matters. A proposition that must be REACHED is under-approximated (INNER), so the
# abstraction never claims a visit that did not happen; one that must be AVOIDED is
# over-approximated (OUTER), so it never misses a violation. R2 is avoided outright and R3 triggers
# an obligation, so both are OUTER. Labelling everything INNER -- as an earlier run did -- is
# optimistic about the avoid-parts, and the grid then CERTIFIES LESS as it refines, because the
# optimism is withdrawn: 50, 35, 11 certified probes at h = 1.0, 0.5, 0.25.
const AP_SEMANTICS =
    Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER, :R2 => MP.OUTER, :R3 => MP.OUTER)
regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
x0 = SVector(-4.0, -7.0)   # initial point a in the paper

println("\nquotient: ", length(quotient.states), " cells; solving their formula ...")
qres = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    regions,
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)
qwin = Set(qres.controllable_set)
println("quotient certifies ", length(qwin), " of ", length(quotient.states), " cells")

# ---------------------------------------------------------
# The grid, on the same domain, answering the same formula
# ---------------------------------------------------------

function grid_cosafe(h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        domain,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    concrete = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]]),
        Dionysos.spot_stepper(φ),
        regions,
        AP_SEMANTICS,
    )
    opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), concrete)
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(h, h))
    MOI.set(opt, MOI.RawOptimizerAttribute("state_grid"), grid)
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(1.0), SVector(1.0)),
    )
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("overapproximation_map"), exact_image_box)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(opt)

    abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
    inner = opt.control_solver.abstract_optimizer
    return (;
        grid,
        mapping = SY.get_state_mapping(abstract_system),
        n_states = SY.get_n_state(abstract_system),
        controllable = Set(inner.controllable_set),
    )
end

function grid_certifies(g, x)
    state = try
        MP.get_state_by_pos(g.mapping, MP.get_pos_by_coord(g.grid, x))
    catch
        return false
    end
    return state !== nothing && state in g.controllable
end

quotient_certifies(x) =
    any(id -> haskey(quotient.states, id) && x in quotient.states[id].set, qwin)
covered(x) = any(s -> x in s.set, values(quotient.states))

box = LazySets.box_approximation(domain)
probes = [
    SVector(x, y) for
    x in range(LazySets.low(box)[1], LazySets.high(box)[1]; length = 45) for
    y in range(LazySets.low(box)[2], LazySets.high(box)[2]; length = 45) if
    SVector(x, y) in domain
]

println(
    "\nGol-Lazar-Belta, their formula, common domain (area ",
    round(LazySets.volume(domain); digits = 2),
    "), ",
    length(probes),
    " probes\n",
)
println(
    rpad("h", 7),
    rpad("grid cells", 12),
    rpad("both", 7),
    rpad("quotient only", 15),
    rpad("GRID ONLY", 11),
    "verdict",
)
for h in [1.0, 0.5, 0.25]
    g = try
        grid_cosafe(h)
    catch err
        println(rpad(h, 7), "FAILED: ", first(split(sprint(showerror, err), "\n")))
        continue
    end
    both, q_only, g_only = 0, 0, 0
    for x in probes
        covered(x) || continue
        q_yes, g_yes = quotient_certifies(x), grid_certifies(g, x)
        q_yes && g_yes ? (both += 1) : q_yes ? (q_only += 1) : g_yes && (g_only += 1)
    end
    println(
        rpad(h, 7),
        rpad(g.n_states, 12),
        rpad(both, 7),
        rpad(q_only, 15),
        rpad(g_only, 11),
        g_only == 0 ? "containment holds" : "*** CONTAINMENT VIOLATED ***",
    )
end
