
module PCLFBisimulationQuotient

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework
const PR = DI.Problem

using JuMP
using LazySets
import HybridSystems
import LinearAlgebra as LA
import Base
import RecipesBase: @recipe, @series
using JLD2

include("geometry_interface.jl")
include("bisimulation_quotient.jl")
include("sublevel_support.jl")

include("cosafe_ltl_problem.jl")

mutable struct OptimizerBisimulationQuotient{T} <: MOI.AbstractOptimizer
    # --- user inputs ---
    bisimulation_quotient_problem::Union{Nothing, PR.BisimulationQuotientProblem}
    pclf::Union{Nothing, PCLF.PCLF}

    level_tol::T
    max_levels::Int
    max_slices::Int

    atol::T
    verbose::Bool

    # --- results ---
    Γ::Union{Nothing, Vector{Float64}}
    D::Union{Nothing, SemiLinearSet}
    bisimulation_quotient::Any
    construction_time_sec::T

    function OptimizerBisimulationQuotient{T}() where {T}
        return new{T}(
            nothing,    # bisimulation_quotient_problem
            nothing,    # pclf
            1e-3,       # level_tol
            200,        # max_levels    
            200,        # max_slices
            1e-3,       # atol
            true,       # verbose
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

function MOI.set(
    model::OptimizerBisimulationQuotient,
    param::MOI.RawOptimizerAttribute,
    value,
)
    name = Symbol(param.name)
    if !hasproperty(model, name)
        error("Unknown optimizer attribute: $(param.name)")
    end
    setproperty!(model, name, value)
    return
end

function MOI.get(model::OptimizerBisimulationQuotient, ::MOI.SolveTimeSec)
    return model.construction_time_sec
end

function MOI.get(model::OptimizerBisimulationQuotient, param::MOI.RawOptimizerAttribute)
    name = Symbol(param.name)
    if !hasproperty(model, name)
        error("Unknown optimizer attribute: $(param.name)")
    end
    return getproperty(model, name)
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

function export_optimizer_jld2(opt::OptimizerBisimulationQuotient, filename::AbstractString)
    jldopen(filename, "w") do f
        f["format_version"] = 1
        return f["optimizer"] = opt
    end
    return nothing
end

function import_optimizer_jld2(filename::AbstractString)
    return jldopen(filename, "r") do f
        v = f["format_version"]
        v == 1 || error("Unsupported optimizer file format_version=$v")
        return f["optimizer"]
    end
end

function export_bisimulation_jld2(
    opt::OptimizerBisimulationQuotient,
    filename::AbstractString,
)
    bisimulation_quotient = opt.bisimulation_quotient
    bisimulation_quotient === nothing && error("No bisimulation quotient computed yet.")

    jldopen(filename, "w") do f
        # versioning for forward compatibility
        f["format_version"] = 1
        f["bisimulation_quotient"] = bisimulation_quotient
        return f
    end
    return nothing
end

function import_bisimulation_jld2(
    filename::AbstractString;
    opt::Union{Nothing, OptimizerBisimulationQuotient} = nothing,
)
    if opt === nothing
        opt = MOI.instantiate(OptimizerBisimulationQuotient)
    end
    jldopen(filename, "r") do f
        v = f["format_version"]
        v == 1 || error("Unsupported bisimulation file format_version=$v")

        bisimulation_quotient = f["bisimulation_quotient"]

        return MOI.set(
            opt,
            MOI.RawOptimizerAttribute("bisimulation_quotient"),
            bisimulation_quotient,
        )
    end
    return opt
end

function MOI.optimize!(opt::OptimizerBisimulationQuotient)
    t_ref = time()

    _validate_model(opt, [:bisimulation_quotient_problem, :pclf])

    prob = opt.bisimulation_quotient_problem

    system = prob.system
    X = prob.region
    regions = prob.observation_regions

    Γ, D = build_levels_and_terminal_set(
        opt.pclf,
        X,
        [_as_hpolytope(R) for R in regions];
        tol = opt.level_tol,
        max_levels = opt.max_levels,
    )
    opt.Γ = Γ
    opt.D = D

    println("Computed levels Γ = ", Γ)
    println("Computed terminal set D = ", D)

    T = bisimulation_pclf(
        system,
        opt.pclf,
        opt.Γ,
        regions;
        verbose = opt.verbose,
        atol = opt.atol,
        max_slices = opt.max_slices,
    )

    opt.bisimulation_quotient = T
    opt.construction_time_sec = time() - t_ref
    return
end

# ============================================================
# Main PCLF bisimulation algorithm
# ============================================================

function bisimulation_pclf(
    f::HybridSystems.HybridSystem,
    pclf::PCLF.PCLF,
    Γ::AbstractVector{<:Real},
    regions;
    verbose::Bool = true,
    atol::Float64 = 0.0,
    max_slices::Union{Nothing, Int} = nothing,
)
    A = extract_mode_matrices(f)

    sublevels = build_sublevel_sequence(pclf, Γ)
    slices = build_slice_sequence(sublevels; atol = atol)

    U = typeof(first(pclf.graph.verts))
    SL = SemiLinearSet
    T = PCBisimulationQuotient{SL, U}(slices)

    initialize_partitions!(T; neutral_obs = 0, terminal_obs = -1)
    refine_partitions_by_observations!(T, regions; terminal_obs = -1, atol = atol)
    initialize_terminal_transitions!(T, pclf)

    incoming_by_dm = group_edges_by_dest_mode(pclf.graph.edges, U)

    N = length(Γ)
    refine_count = 0
    for i in 1:min(max_slices, N)
        verbose && println("Current slice = $i")

        for ((d, m), source_nodes) in incoming_by_dm
            verbose && println("  destination/mode = ($d, $m)")

            target_ids = [
                qid for qid in get(T.part_ids, d, Int[]) if
                haskey(T.states, qid) && T.states[qid].slice == i
            ]

            isempty(target_ids) && continue

            for qid in target_ids
                haskey(T.states, qid) || continue
                q = T.states[qid]

                pre_parts = preimage_linear_parts(q.set, A[m])
                isempty(pre_parts) && continue

                for s in source_nodes
                    source_ids = [
                        pid for pid in get(T.part_ids, s, Int[]) if
                        haskey(T.states, pid) && T.states[pid].slice > i
                    ]

                    for pid in source_ids
                        haskey(T.states, pid) || continue
                        qsrc = T.states[pid]

                        if any(t -> t[1] == m, qsrc.next) # to investigate
                            continue
                        end

                        refined = refine_one_state!(T, pid, pre_parts, m, qid; atol = atol)
                        refined && (refine_count += 1)

                        verbose &&
                            refined &&
                            refine_count % 50 == 0 &&
                            println("    refinement: $refine_count")
                    end
                end
            end
        end
    end

    return T
end

# ============================================================
# Switched system helpers
# ============================================================

# could be move to system
function extract_mode_matrices(f::HybridSystems.HybridSystem)
    RMs = f.resetmaps
    A = Vector{Matrix{Float64}}(undef, length(RMs))
    for (i, rm) in enumerate(RMs)
        if isa(rm, AbstractMatrix)
            A[i] = Array(rm)
        elseif :A in fieldnames(typeof(rm))
            A[i] = Array(getfield(rm, :A))
        else
            error("Cannot extract matrix from resetmap of type $(typeof(rm)).")
        end
    end
    return A
end

# ============================================================
# Slice generation
# ============================================================

function build_sublevel_sequence(pclf::PCLF.PCLF, Γ::AbstractVector{<:Real})
    U = typeof(first(pclf.graph.verts))
    sublevels = Dict{U, Vector{Poly}}()
    for s in pclf.graph.verts
        piece = pclf.pieces[s]
        sublevels[s] = [_as_hpolytope(PCLF.get_sublevel_set(piece, Float64(γ))) for γ in Γ]
    end
    return sublevels
end

function build_slice_sequence(
    sublevels::Dict{U, Vector{Poly}};
    atol::Float64 = 0.0,
) where {U}
    slices = Dict{U, Vector{SemiLinearSet}}()

    for (s, Ps) in sublevels
        Ns = length(Ps)
        local_slices = Vector{SemiLinearSet}(undef, Ns)

        local_slices[1] = SemiLinearSet(Ps[1])

        for i in 2:Ns
            diff_parts = set_difference_decompose(Ps[i], Ps[i - 1]; atol = atol)
            local_slices[i] = SemiLinearSet(diff_parts)
        end
        slices[s] = local_slices
    end
    return slices
end

# ============================================================
# Initialize partitions and transitions
# ============================================================

function initialize_partitions!(
    T::PCBisimulationQuotient{SemiLinearSet, U};
    neutral_obs::Int = 0,
    terminal_obs::Int = -1,
) where {U}
    for (s, slice_list) in T.slices
        for (i, Sset) in enumerate(slice_list)
            obs = (i == 1) ? terminal_obs : neutral_obs
            if is_nonempty_set(Sset)
                add_state!(T, s, Sset, obs, i)
            end
        end
    end
    return T
end

function refine_state_by_observation!(
    T::PCBisimulationQuotient{SemiLinearSet, U},
    qid::Int,
    R::Poly,
    new_obs::Int;
    terminal_obs::Int = -1,
    atol::Float64 = 0.0,
) where {U}
    haskey(T.states, qid) || return false
    q = T.states[qid]

    q.obs == terminal_obs && return false

    old_next = copy(q.next)
    old_node = q.node
    old_slice = q.slice
    old_obs = q.obs

    inside_parts = Poly[]
    outside_parts = Poly[]
    touched = false

    for Q in q.set
        I = set_intersection(Q, R)
        if is_nonempty_set(I)
            touched = true
            push!(inside_parts, I)

            D = set_difference_decompose(Q, R; atol = atol)
            if !isempty(D)
                append!(outside_parts, D)
            end
        else
            push!(outside_parts, Q)
        end
    end

    touched || return false

    remove_state!(T, qid)

    if !isempty(outside_parts)
        out_id = add_state!(T, old_node, SemiLinearSet(outside_parts), old_obs, old_slice)
        T.states[out_id].next = old_next
    end

    if !isempty(inside_parts)
        in_id = add_state!(T, old_node, SemiLinearSet(inside_parts), new_obs, old_slice)
        T.states[in_id].next = old_next
    end

    return true
end

function refine_partitions_by_observations!(
    T::PCBisimulationQuotient{SemiLinearSet, U},
    regions;
    terminal_obs::Int = -1,
    atol::Float64 = 0.0,
) where {U}
    Rh = [_as_hpolytope(R) for R in regions]

    for (obs, R) in enumerate(Rh)
        qids = copy(collect(keys(T.states)))
        for qid in qids
            haskey(T.states, qid) || continue
            refine_state_by_observation!(
                T,
                qid,
                R,
                obs;
                terminal_obs = terminal_obs,
                atol = atol,
            )
        end
    end
    return T
end

function initialize_terminal_transitions!(T::PCBisimulationQuotient, pclf::PCLF.PCLF)
    terminal_by_node = Dict{Any, Vector{Int}}()

    for (qid, q) in T.states
        if q.slice == 1
            get!(terminal_by_node, q.node, Int[])
            push!(terminal_by_node[q.node], qid)
        end
    end

    for (s, d, m) in pclf.graph.edges
        if haskey(terminal_by_node, s) && haskey(terminal_by_node, d)
            for qs in terminal_by_node[s], qd in terminal_by_node[d]
                add_transition!(T, qs, Int(m), qd)
            end
        end
    end

    return T
end

# ============================================================
# Refinement
# ============================================================

function refine_one_state!(
    T::PCBisimulationQuotient{SemiLinearSet, U},
    qid::Int,
    pre_parts::AbstractVector,
    mode::Int,
    target_qid::Int;
    atol::Float64 = 0.0,
) where {U}
    haskey(T.states, qid) || return false
    q = T.states[qid]

    old_next = copy(q.next)
    old_obs = q.obs
    old_node = q.node
    old_slice = q.slice

    inside_parts = Poly[]
    outside_parts = Poly[]
    touched = false

    for Q in q.set
        Qrem = SemiLinearSet(Q)

        for preP in pre_parts
            I = set_intersection(Qrem, preP)
            if !is_nonempty_set(I)
                continue
            end

            touched = true
            append!(inside_parts, I)

            Qrem = set_difference_decompose(Qrem, preP; atol = atol)
            if !is_nonempty_set(Qrem)
                break
            end
        end

        if is_nonempty_set(Qrem)
            append!(outside_parts, Qrem)
        end
    end

    touched || return false

    remove_state!(T, qid)

    if !isempty(outside_parts)
        diff_id = add_state!(T, old_node, SemiLinearSet(outside_parts), old_obs, old_slice)
        T.states[diff_id].next = copy(old_next)
    end

    inter_id = add_state!(T, old_node, SemiLinearSet(inside_parts), old_obs, old_slice)
    T.states[inter_id].next = copy(old_next)
    add_transition!(T, inter_id, mode, target_qid)

    return true
end

# ============================================================
# Helpers
# ============================================================

function group_edges_by_dest_mode(edges, ::Type{U}) where {U}
    grouped = Dict{Tuple{U, Int}, Vector{U}}()
    for (s, d, m) in edges
        push!(get!(grouped, (d, Int(m)), U[]), s)
    end
    return grouped
end

end # module
