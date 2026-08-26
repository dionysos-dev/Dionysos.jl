
module PCLFBisimulationQuotient

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework
const ST = DI.System
const PR = DI.Problem
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems

using JuMP

import LazySets
import Polyhedra

import HybridSystems
import LinearAlgebra as LA

import Base
import RecipesBase: @recipe, @series

const Poly = LazySets.HPolytope

include("bisimulation_quotient.jl")
include("sublevel_support.jl")

include("cosafe_ltl_problem.jl")

"""
    OptimizerBisimulationQuotient{T} <: Dionysos.Optim.AbstractDionysosOptimizer

Builds the path-complete-Lyapunov bisimulation quotient of a switched system.

Set `"bisimulation_quotient_problem"` and `"pclf"`; read back `"bisimulation_quotient"`,
`"D"` (the terminal set) and `"construction_time_sec"`.

`atol` is the inset used when one polytope is cut out of another. It is not a mere
tolerance: every cut moves the cutting plane inward by `atol`, so the quotient loses a thin
shell at each one and the loss compounds over the thousands of cuts a refinement sweep
performs. Measured on a 7-level, 2-mode example, the fraction of the domain left
uncovered tracks `atol` linearly — 3.7% at `1e-3`, 0.39% at `1e-4`, 0.004% at `1e-6` —
while the cost of tightening it is slight (a fifth more time from `1e-3` to `1e-6`, and a
few percent more states). It cannot go to zero: cutting exactly produces degenerate,
flat pieces that survive the emptiness test and make the construction diverge.
"""
mutable struct OptimizerBisimulationQuotient{T} <: OP.AbstractDionysosOptimizer
    # --- user inputs ---
    bisimulation_quotient_problem::Union{Nothing, PR.BisimulationQuotientProblem}
    pclf::Union{Nothing, PCLF.PCLF}

    level_tol::T
    max_levels::Int
    ΓX::Union{Nothing, Float64}
    nb_levels::Union{Nothing, Int}
    max_slices::Int
    polyhedra_backend::Any

    atol::T
    print_level::Int

    # --- results ---
    Γ::Union{Nothing, Vector{Float64}}
    D::Union{Nothing, UT.SemiLinearSet}
    bisimulation_quotient::Any
    construction_time_sec::T

    function OptimizerBisimulationQuotient{T}() where {T}
        return new{T}(
            nothing,    # bisimulation_quotient_problem
            nothing,    # pclf
            1e-2,       # level_tol
            200,        # max_levels
            nothing,    # ΓX
            nothing,    # nb_levels
            200,        # max_slices
            nothing,    # polyhedra_backend
            1e-6,       # atol
            1,          # print_level
            nothing,    # Γ
            nothing,    # D
            nothing,    # bisimulation_quotient
            zero(T),    # construction_time_sec
        )
    end
end

OptimizerBisimulationQuotient() = OptimizerBisimulationQuotient{Float64}()

MOI.is_empty(opt::OptimizerBisimulationQuotient) =
    opt.bisimulation_quotient_problem === nothing

function MOI.get(model::OptimizerBisimulationQuotient, ::MOI.SolveTimeSec)
    return model.construction_time_sec
end

function reset!(model::OptimizerBisimulationQuotient)
    model.Γ = nothing
    model.D = nothing
    model.bisimulation_quotient = nothing
    model.construction_time_sec = 0.0
    return model
end

function _validate_model(
    model::OptimizerBisimulationQuotient,
    required_fields::Vector{Symbol},
)
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set `$(field)`. Missing required field in OptimizerBisimulationQuotient.",
            )
        end
    end
end

function MOI.optimize!(opt::OptimizerBisimulationQuotient)
    t_ref = time()

    _validate_model(opt, [:bisimulation_quotient_problem, :pclf])

    prob = opt.bisimulation_quotient_problem

    system = prob.system
    # the level builder works on HPolytopes; accept any polyhedral region
    # (e.g. a Hyperrectangle) by converting at the boundary
    X = UT._as_hpolytope(prob.state_set)
    regions = prob.observation_regions

    Γ, D = build_levels_and_terminal_set(
        opt.pclf,
        X,
        [UT._as_hpolytope(R) for R in regions];
        tol = opt.level_tol,
        max_levels = opt.max_levels,
        ΓX = opt.ΓX,
        nb_levels = opt.nb_levels,
    )
    opt.Γ = Γ
    opt.D = D

    if opt.print_level >= 1
        println("Computed levels Γ = ", Γ)
        println("Computed terminal set D = ", D)
    end

    quotient = bisimulation_pclf(
        system,
        opt.pclf,
        opt.Γ,
        regions;
        verbose = opt.print_level >= 1,
        atol = opt.atol,
        max_slices = opt.max_slices,
    )

    opt.bisimulation_quotient = quotient
    opt.construction_time_sec = time() - t_ref
    return
end

# ============================================================
# Main PCLF bisimulation algorithm
# ============================================================

"""
    refine_sources!(quotient, source_nodes, target_qid, pre_parts, mode, dest_node, slice; atol)

Split every state of `source_nodes` that partly reaches `target_qid` under `mode`, and return
how many were split.

Sources are read from the quotient inside the loop rather than collected up front: splitting a
state removes it and puts its pieces back, and a piece may still have to be split against a
later target. A source that already reaches `dest_node` under `mode` is left alone -- it has
been split for that destination already, and splitting it again would only fragment the
partition without adding transitions.
"""
function refine_sources!(
    quotient::PCBisimulationQuotient,
    source_nodes,
    target_qid::Int,
    pre_parts::AbstractVector,
    mode::Int,
    dest_node,
    slice::Int;
    atol::Float64 = 1e-6,
)
    refined = 0
    for s in source_nodes
        source_ids = [
            pid for pid in get(quotient.part_ids, s, Int[]) if
            haskey(quotient.states, pid) && quotient.states[pid].slice > slice
        ]
        for pid in source_ids
            haskey(quotient.states, pid) || continue
            any(
                t -> t[1] == mode && quotient.states[t[2]].node == dest_node,
                quotient.states[pid].next,
            ) && continue
            refine_one_state!(quotient, pid, pre_parts, mode, target_qid; atol = atol) &&
                (refined += 1)
        end
    end
    return refined
end

"""
    bisimulation_pclf(system, pclf, Γ, regions; verbose, atol, max_slices)

Build the path-complete-Lyapunov bisimulation quotient of a discrete switched `system`.

The construction works inward from the level sets of `pclf`. Each Lyapunov piece is cut at
every level of `Γ` and the nested sublevel sets are turned into disjoint shells, one
abstract state per shell; the innermost shell is the terminal set. States are then split
so that no one of them straddles an observation region, and finally refined slice by slice:
a source state is split against the preimage of a target so that the part that reaches the
target -- and only that part -- carries the transition. Working outward one slice at a time
is what makes the result a bisimulation rather than an over-approximation, since a state is
only ever split against targets whose own refinement is already settled.

`atol` is the inset applied whenever one polytope is cut out of another; see
[`OptimizerBisimulationQuotient`](@ref) for what it costs. It must be strictly positive --
cutting exactly leaves degenerate pieces that survive the emptiness test and make the
refinement diverge.
"""
function bisimulation_pclf(
    system::HybridSystems.HybridSystem,
    pclf::PCLF.PCLF,
    Γ::AbstractVector{<:Real},
    regions;
    verbose::Bool = true,
    atol::Float64 = 1e-6,
    max_slices::Union{Nothing, Int} = nothing,
)
    A = ST.mode_matrices(system)

    sublevels = build_sublevel_sequence(pclf, Γ)
    slices = build_slice_sequence(sublevels; atol = atol)

    U = typeof(first(pclf.graph.verts))
    SL = UT.SemiLinearSet
    quotient = PCBisimulationQuotient{SL, U}(slices)

    initialize_partitions!(quotient; neutral_obs = 0, terminal_obs = -1)
    refine_partitions_by_observations!(quotient, regions; terminal_obs = -1, atol = atol)
    initialize_terminal_transitions!(quotient, pclf)

    sources_by_dest_mode = group_edges_by_dest_mode(pclf.graph.edges, U)

    N = length(Γ)
    refine_count = 0
    for i in 1:min(max_slices, N)
        verbose && println("Current slice = $i")

        for ((d, m), source_nodes) in sources_by_dest_mode
            verbose && println("  destination/mode = ($d, $m)")

            target_ids = [
                qid for qid in get(quotient.part_ids, d, Int[]) if
                haskey(quotient.states, qid) && quotient.states[qid].slice == i
            ]

            isempty(target_ids) && continue

            for qid in target_ids
                haskey(quotient.states, qid) || continue

                pre_parts = UT.preimage_linear_parts(quotient.states[qid].set, A[m])
                isempty(pre_parts) && continue

                before = refine_count
                refine_count += refine_sources!(
                    quotient,
                    source_nodes,
                    qid,
                    pre_parts,
                    m,
                    d,
                    i;
                    atol = atol,
                )
                verbose &&
                    refine_count > before &&
                    div(before, 50) != div(refine_count, 50) &&
                    println("    refinement: $refine_count")
            end
        end
    end

    return quotient
end

# ============================================================
# Switched system helpers
# ============================================================

# ============================================================
# Slice generation
# ============================================================

"""
    build_sublevel_sequence(pclf, Γ)

The sublevel sets of every Lyapunov piece, one per level of `Γ`, keyed by graph node.

A piece may be a single polytope or a union of them, so the values are heterogeneous; the
caller turns them into shells.
"""
function build_sublevel_sequence(pclf::PCLF.PCLF, Γ::AbstractVector{<:Real})
    U = typeof(first(pclf.graph.verts))
    sublevels = Dict{U, Vector{Union{Poly, UT.SemiLinearSet}}}()
    for s in pclf.graph.verts
        piece = pclf.pieces[s]
        sublevels[s] = [PCLF.get_sublevel_set(piece, Float64(γ)) for γ in Γ]
    end
    return sublevels
end

"""
    build_slice_sequence(sublevels; atol)

Turn each node's nested sublevel sets into the disjoint shells the abstract states are built
from: the innermost sublevel set, then the difference of each level with the one inside it.

The shells of a node cover its outermost sublevel set up to the `atol` inset each cut leaves
behind, and they are pairwise disjoint -- which is what lets a volume over the quotient be a
plain sum within a node.
"""
function build_slice_sequence(
    sublevels::Dict{U, Vector{Union{Poly, UT.SemiLinearSet}}};
    atol::Float64 = 1e-6,
) where {U}
    slices = Dict{U, Vector{UT.SemiLinearSet}}()

    for (s, Ps) in sublevels
        Ns = length(Ps)
        local_slices = Vector{UT.SemiLinearSet}(undef, Ns)

        if Ps[1] isa Poly
            local_slices[1] = UT.semilinear_set([Ps[1]])
        else
            local_slices[1] = Ps[1]
        end

        for i in 2:Ns
            diff_parts = UT.set_difference_decompose(Ps[i], Ps[i - 1]; atol = atol)
            local_slices[i] = UT.semilinear_set(diff_parts)
        end
        slices[s] = local_slices
    end
    return slices
end

# ============================================================
# Initialize partitions and transitions
# ============================================================

"""
    initialize_partitions!(quotient; neutral_obs, terminal_obs)

Seed one abstract state per non-empty shell.

The innermost shell of every node is the terminal set and is marked `terminal_obs`; everything
else starts unobserved and is given its real observation by
[`refine_partitions_by_observations!`](@ref).
"""
function initialize_partitions!(
    quotient::PCBisimulationQuotient{UT.SemiLinearSet, U};
    neutral_obs::Int = 0,
    terminal_obs::Int = -1,
) where {U}
    for (s, slice_list) in quotient.slices
        for (i, Sset) in enumerate(slice_list)
            obs = (i == 1) ? terminal_obs : neutral_obs
            if !isempty(Sset)
                add_state!(quotient, s, Sset, obs, i)
            end
        end
    end
    return quotient
end

"""
    refine_state_by_observation!(quotient, qid, R, new_obs; terminal_obs, atol)

Split state `qid` along region `R`, giving the part inside `R` the observation `new_obs` and
leaving the rest as it was. Returns whether anything was split.

Both halves inherit the outgoing transitions of the original: a split refines *where* the
state is, not what it can do. Terminal states are left alone -- their observation is what
makes them terminal.
"""
function refine_state_by_observation!(
    quotient::PCBisimulationQuotient{UT.SemiLinearSet, U},
    qid::Int,
    R::Poly,
    new_obs::Int;
    terminal_obs::Int = -1,
    atol::Float64 = 1e-6,
) where {U}
    haskey(quotient.states, qid) || return false
    q = quotient.states[qid]

    q.obs == terminal_obs && return false

    old_next = copy(q.next)
    old_node = q.node
    old_slice = q.slice
    old_obs = q.obs

    inside_parts = Poly[]
    outside_parts = Poly[]
    touched = false

    for Q in q.set.array
        # `poly_intersection` prunes: an empty intersection is an `EmptySet`
        I = UT.poly_intersection(Q, R)
        if !(I isa LazySets.EmptySet)
            touched = true
            push!(inside_parts, I)

            D = UT.set_difference_decompose(Q, R; atol = atol) # maybe something to change here??
            if !isempty(D)
                append!(outside_parts, D)
            end
        else
            push!(outside_parts, Q)
        end
    end

    touched || return false

    remove_state!(quotient, qid)

    if !isempty(outside_parts)
        out_id = add_state!(
            quotient,
            old_node,
            UT.semilinear_set(outside_parts),
            old_obs,
            old_slice,
        )
        quotient.states[out_id].next = old_next
    end

    if !isempty(inside_parts)
        in_id = add_state!(
            quotient,
            old_node,
            UT.semilinear_set(inside_parts),
            new_obs,
            old_slice,
        )
        quotient.states[in_id].next = old_next
    end

    return true
end

"""
    refine_partitions_by_observations!(quotient, regions; terminal_obs, atol)

Split states until none of them straddles the boundary of an observation region.

Every state must have a single observation for the quotient to respect the labelling a
specification is written against, so each region is applied in turn to every state alive at
that point -- including the states produced by the previous region.
"""
function refine_partitions_by_observations!(
    quotient::PCBisimulationQuotient{UT.SemiLinearSet, U},
    regions;
    terminal_obs::Int = -1,
    atol::Float64 = 1e-6,
) where {U}
    region_polytopes = [UT._as_hpolytope(R) for R in regions]

    for (obs, R) in enumerate(region_polytopes)
        qids = copy(collect(keys(quotient.states)))
        for qid in qids
            haskey(quotient.states, qid) || continue
            refine_state_by_observation!(
                quotient,
                qid,
                R,
                obs;
                terminal_obs = terminal_obs,
                atol = atol,
            )
        end
    end
    return quotient
end

"""
    initialize_terminal_transitions!(quotient, pclf)

Connect the terminal states along the edges of the Lyapunov graph.

The innermost shell is where the certificate stops distinguishing states, so any terminal
state is taken to reach any terminal state of a node the graph has an edge to. These are the
only transitions not discovered by refinement.
"""
function initialize_terminal_transitions!(quotient::PCBisimulationQuotient, pclf::PCLF.PCLF)
    terminal_by_node = Dict{Any, Vector{Int}}()

    for (qid, q) in quotient.states
        if q.slice == 1
            get!(terminal_by_node, q.node, Int[])
            push!(terminal_by_node[q.node], qid)
        end
    end

    for (s, d, m) in pclf.graph.edges
        if haskey(terminal_by_node, s) && haskey(terminal_by_node, d)
            for qs in terminal_by_node[s], qd in terminal_by_node[d]
                add_transition!(quotient, qs, Int(m), qd)
            end
        end
    end

    return quotient
end

# ============================================================
# Refinement
# ============================================================

"""
    refine_one_state!(quotient, qid, pre_parts, mode, target_qid; atol)

Split source state `qid` against `pre_parts`, the preimage of a target under `mode`, so that
the part which lands in the target carries the transition and the rest does not. Returns
whether anything was split.

This is the step that makes the quotient exact: without it a state that only partly reaches
the target would carry the transition wholesale, and the abstraction would admit behaviour the
system does not have. Both halves inherit the original outgoing transitions; only the part
inside the preimage gains the new one.
"""
function refine_one_state!(
    quotient::PCBisimulationQuotient{UT.SemiLinearSet, U},
    qid::Int,
    pre_parts::AbstractVector,
    mode::Int,
    target_qid::Int;
    atol::Float64 = 1e-6,
) where {U}
    haskey(quotient.states, qid) || return false
    q = quotient.states[qid]

    old_next = copy(q.next)
    old_obs = q.obs
    old_node = q.node
    old_slice = q.slice

    inside_parts = Poly[]
    outside_parts = Poly[]
    touched = false

    for Q in q.set.array
        remainder = UT.semilinear_set(Q)

        for preP in pre_parts
            hits = UT.poly_intersection_parts(remainder, preP)
            if isempty(hits)
                continue
            end

            touched = true
            append!(inside_parts, hits)

            # decompose only keeps nonempty pieces, so structural emptiness
            # of the remainder is exact
            parts = UT.set_difference_decompose(remainder, preP; atol = atol)
            remainder = UT.semilinear_set(parts)
            if isempty(remainder.array)
                break
            end
        end

        append!(outside_parts, remainder.array)
    end

    touched || return false

    remove_state!(quotient, qid)

    if !isempty(outside_parts)
        diff_id = add_state!(
            quotient,
            old_node,
            UT.semilinear_set(outside_parts),
            old_obs,
            old_slice,
        )
        quotient.states[diff_id].next = copy(old_next)
    end

    inter_id =
        add_state!(quotient, old_node, UT.semilinear_set(inside_parts), old_obs, old_slice)
    quotient.states[inter_id].next = copy(old_next)
    add_transition!(quotient, inter_id, mode, target_qid)

    return true
end

# ============================================================
# Helpers
# ============================================================

"""
    group_edges_by_dest_mode(edges, U)

Index the Lyapunov graph edges by `(destination, mode)`, giving the source nodes for each.

Refinement walks targets and asks which nodes can reach them, which is the reverse of how the
edges are stored.
"""
function group_edges_by_dest_mode(edges, ::Type{U}) where {U}
    grouped = Dict{Tuple{U, Int}, Vector{U}}()
    for (s, d, m) in edges
        push!(get!(grouped, (d, Int(m)), U[]), s)
    end
    return grouped
end

end # module
