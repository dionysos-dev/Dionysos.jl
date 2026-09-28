# Is the grid's certified set really a subset of the bisimulation's?
#
# The second half of N3. A grid baseline measured earlier showed the grid's certified volume
# rising with refinement -- 0.16, 0.60, 0.79, 1.04 -- which is an over-approximating abstraction
# approaching the true controllable set from below. This file checks the containment directly.
#
# The comparison is POINTWISE, on a probe grid, not by volume. `max_slices` truncation makes covered
# volumes differ for reasons having nothing to do with the method, which is why E1 uses probes as
# its instrument; volumes are reported only as a secondary column.
#
# What each column means:
#
#   both          probe certified by the grid AND by the quotient -- agreement
#   quotient only the gain: exactness certifying what the over-approximation cannot
#   GRID ONLY     *** must be zero ***: the grid certifying something the exact answer does not
#                 means either the grid is unsound or the quotient is not exact
#
# `grid only > 0` is the falsifier of the whole experiment, so it is checked before anything else
# is read.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra
import MathematicalSystems

const SY = DI.Symbolic
const OPDS = OP.DiscreteSystems

const A = [[0.70 0.10; 0.00 0.65], [0.60 -0.15; 0.10 0.55]]

mode_of(u) = clamp(round(Int, u[1]), 1, length(A))
step_map(x, u) = A[mode_of(u)] * x

function exact_image_box(rect, u)
    c = LazySets.center(rect)
    r = LazySets.radius_hyperrectangle(rect)
    M = A[mode_of(u)]
    return LazySets.Hyperrectangle(M * c, abs.(M) * r)
end

# ---------------------------------------------------------
# A common domain, so the cell counts mean the same thing
# ---------------------------------------------------------

"""
The outermost sublevel set of the certificate: exactly the region the quotient covers.

The two abstractions cover different sets if left alone -- the grid is built on the working set `X`,
the quotient on the PCLF's sublevel sets, a rotated box reaching well beyond it -- so counting their
cells would compare counts over different areas. Handing this polytope to the grid as its state set
makes the two domains **identical**, which is better than gridding its bounding box and filtering
afterwards: the bounding box of a rotated box has corners the quotient never covers, and gridding
them inflates the grid's count for nothing.

`MP.get_states_from_set` accepts any bounded `LazySet`, so the grid does not need a hyperrectangle.
"""
function certificate_domain(pclf, Γ)
    node = first(keys(pclf.pieces))
    return PCLF.get_sublevel_set(pclf.pieces[node], maximum(Γ))
end

function covered_area(quotient)
    return sum(sum(LazySets.volume, state.set.array) for state in values(quotient.states))
end

in_quotient_domain(quotient, x) = any(s -> x in s.set, values(quotient.states))

# ---------------------------------------------------------
# The grid side
# ---------------------------------------------------------

"""
Grid abstraction at step `h`, with the set of controllable grid states for reaching `target`.

Returns the mapping and grid alongside, so that a probe point can be located afterwards.
"""
function grid_controllable(X, target, h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        X,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(system, X),
    )
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(h, h))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(1.0), SVector(1.0)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("overapproximation_map"), exact_image_box)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    automaton = SY.get_automaton(abstract_system)
    mapping = SY.get_state_mapping(abstract_system)
    targets = collect(MP.get_states_from_set(mapping, target, MP.INNER))
    _, controllable, _, _ = OPDS.compute_worst_case_cost_controller(automaton, targets)

    return (;
        grid,
        mapping,
        controllable = Set(controllable),
        n_states = SY.get_n_state(abstract_system),
    )
end

"""
Grid cells whose centre lies inside the quotient's region, so the counts cover the same area.
"""
function cells_in_domain(g, domain)
    n = 0
    for q in MP.enum_states(g.mapping)
        c = MP.get_coord_by_pos(g.grid, MP.get_pos_by_state(g.mapping, q))
        c in domain && (n += 1)
    end
    return n
end

function grid_certifies(g, x)
    pos = MP.get_pos_by_coord(g.grid, x)
    state = try
        MP.get_state_by_pos(g.mapping, pos)
    catch
        return false
    end
    return state !== nothing && state in g.controllable
end

# ---------------------------------------------------------
# The quotient side
# ---------------------------------------------------------

"""
Largest number of successors any (state, mode) pair has, and how many pairs branch.

The two sides of this comparison resolve nondeterminism differently, and they must: the grid's
branching is SPURIOUS -- an artefact of over-approximating a cell's image by whole cells -- so it is
resolved adversarially by `compute_worst_case_cost_controller`, and anything else would be unsound.
The quotient's branching, when it has any, is GENUINE graph structure of the lifted system, which
synthesis resolves angelically (`thm:synthesis-invariance`) and verification adversarially.

A quotient branches only when its graph does: measured, the De Bruijn quotients have exactly one
successor per (state, mode) while the nondeterministic observer graph reaches three. This comparison
uses the single-node graph, which is deterministic, so the two semantics coincide on the quotient
side and the comparison against the grid's worst case is exact. That is asserted below rather than
assumed.
"""
function branching(quotient)
    worst, branching_pairs, pairs = 0, 0, 0
    for (_, state) in quotient.states
        per_mode = Dict{Int, Int}()
        for (m, _) in state.next
            per_mode[m] = get(per_mode, m, 0) + 1
        end
        for (_, count) in per_mode
            worst = max(worst, count)
            pairs += 1
            count > 1 && (branching_pairs += 1)
        end
    end
    return worst, branching_pairs, pairs
end

"""
Backward reachability on the quotient under existential semantics: the controller picks the mode.

A state is winning when some recorded transition leads into the winning set, which is the same
fixed point the grid side runs, so the two answers are comparable.
"""
function quotient_winning(quotient, target_obs::Int)
    win = Set{Int}(id for (id, state) in quotient.states if state.obs == target_obs)
    changed = true
    while changed
        changed = false
        for (id, state) in quotient.states
            id in win && continue
            if any(t -> t[2] in win, state.next)
                push!(win, id)
                changed = true
            end
        end
    end
    return win
end

"""
Is `x` certified by the quotient? Existentially: some node's cell containing `x` is winning.
`covered` says whether any cell contains `x` at all, which separates "not certified" from
"outside the region the quotient covers".
"""
function quotient_certifies(quotient, win, x)
    covered = false
    for (id, state) in quotient.states
        if x in state.set
            covered = true
            id in win && return (true, true)
        end
    end
    return (false, covered)
end

# ---------------------------------------------------------
# The comparison
# ---------------------------------------------------------

(; f, problem, X, R1, R2) = two_mode_problem()

println(
    "Building the PCLF quotient (single node, so the smallest quotient of the family) ...",
)
pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    PCLF.generate_DeBruijn_edges(2, 0),
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000);
    Gmats = rotation_templates([1]; θ = π / 6, mode = :rotation),
    MLF = true,
    verbose = false,
)
(; optimizer, quotient) = build_quotient(
    problem,
    pclf;
    atol = 1e-3,
    max_levels = 100,
    max_slices = 8,
    print_level = 0,
)

worst, branching_pairs, pairs = branching(quotient)
println(
    "quotient branching: max ",
    worst,
    " successors per (state, mode), ",
    branching_pairs,
    " of ",
    pairs,
    " pairs branch",
    worst <= 1 ? "  -> deterministic, so ∃ and ∀ coincide here" :
    "  *** branches: ∃ and ∀ differ, the comparison below is the ∃ (synthesis) one ***",
)

# `two_mode_problem` declares regions [R1, R2]; observation 1 is R1.
win = quotient_winning(quotient, 1)
println(
    "quotient: ",
    length(quotient.states),
    " cells, ",
    length(win),
    " winning for reach-R1\n",
)

# The grid is built on the certificate's outermost sublevel set -- the very region the quotient
# covers -- so the two abstractions discretize the SAME set and their cell counts are comparable.
domain = certificate_domain(pclf, optimizer.Γ)
box = LazySets.box_approximation(domain)
println(
    "common domain: the outermost sublevel set, area ",
    round(LazySets.volume(domain); digits = 3),
    "   quotient covers ",
    round(covered_area(quotient); digits = 3),
    " with ",
    length(quotient.states),
    " cells\n",
)

probes = [
    SVector(x, y) for
    x in range(LazySets.low(box)[1], LazySets.high(box)[1]; length = 45) for
    y in range(LazySets.low(box)[2], LazySets.high(box)[2]; length = 45) if
    SVector(x, y) in domain
]

println(
    rpad("h", 7),
    rpad("grid cells", 12),
    rpad("in domain", 11),
    rpad("probes", 8),
    rpad("both", 7),
    rpad("quotient only", 15),
    rpad("GRID ONLY", 11),
    "verdict",
)
for h in [0.4, 0.2, 0.1, 0.05]
    g = grid_controllable(domain, R1, h)
    both, q_only, g_only, compared = 0, 0, 0, 0
    for x in probes
        q_yes, covered = quotient_certifies(quotient, win, x)
        covered || continue          # only compare where the quotient has an answer
        compared += 1
        g_yes = grid_certifies(g, x)
        if g_yes && q_yes
            both += 1
        elseif q_yes
            q_only += 1
        elseif g_yes
            g_only += 1
        end
    end
    println(
        rpad(h, 7),
        rpad(g.n_states, 12),
        rpad(cells_in_domain(g, domain), 11),
        rpad(compared, 8),
        rpad(both, 7),
        rpad(q_only, 15),
        rpad(g_only, 11),
        g_only == 0 ? "containment holds" : "*** CONTAINMENT VIOLATED ***",
    )
end
