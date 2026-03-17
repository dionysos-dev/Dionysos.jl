
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

    obs_partition::Any
    Γ::Union{Nothing, Vector{Float64}}
    num_levels::Union{Nothing, Int}
    tau::Union{Nothing, Float64}
    max_slices::Union{Nothing, Int} # debug

    atol::T
    verbose::Bool

    # --- results ---
    bisimulation_quotient::Any
    abstraction_construction_time_sec::T

    function OptimizerBisimulationQuotient{T}() where {T}
        return new{T}(
            nothing,    # bisimulation_quotient_problem
            nothing,    # pclf
            nothing,    # obs_partition
            nothing,    # Γ
            nothing,    # num_levels
            nothing,    # tau
            nothing,    # max_slices
            zero(T),    # atol
            true,       # verbose
            nothing,    # bisimulation_quotient
            zero(T),    # solve time
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
    return model.abstraction_construction_time_sec
end

function MOI.get(model::OptimizerBisimulationQuotient, param::MOI.RawOptimizerAttribute)
    name = Symbol(param.name)
    if !hasproperty(model, name)
        error("Unknown optimizer attribute: $(param.name)")
    end
    return getproperty(model, name)
end

function reset!(model::OptimizerBisimulationQuotient)
    model.bisimulation_quotient = nothing
    model.abstraction_construction_time_sec = 0.0
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
    D = prob.terminal_region
    regions = prob.observation_regions

    opt.obs_partition = build_observation_partition(
        _as_hpolytope(X),
        _as_hpolytope(D),
        [_as_hpolytope(R) for R in regions];
        atol = opt.atol,
    )

    # Automatic Γ if user did not provide it
    if isnothing(opt.Γ)
        X isa Hyperrectangle || error(
            "Automatic Γ construction currently expects `region` to be a Hyperrectangle.",
        )
        D isa Hyperrectangle || error(
            "Automatic Γ construction currently expects `terminal_region` to be a Hyperrectangle.",
        )

        Γ_auto, τ_auto, ΓD, ΓX =
            build_levels_from_problem(opt.pclf, X, D; num_levels = opt.num_levels)

        opt.Γ = Γ_auto
        opt.tau = τ_auto

        opt.verbose && println("Automatic levels selected:")
        opt.verbose && println("  ΓD = $ΓD")
        opt.verbose && println("  ΓX = $ΓX")
        opt.verbose && println("  τ  = $τ_auto")
        opt.verbose && println("  Γ  = $(opt.Γ)")
    end

    T = bisimulation_pclf(
        system,
        opt.pclf,
        opt.Γ,
        opt.obs_partition;
        verbose = opt.verbose,
        atol = opt.atol,
        max_slices = opt.max_slices,
    )

    opt.bisimulation_quotient = T
    opt.abstraction_construction_time_sec = time() - t_ref
    return
end

# ============================================================
# Observation partition
# ============================================================

"""
    build_observation_partition(X, D, regions; neutral_obs = 0, terminal_obs = -1, atol = 0.0)

Build
    P_X = { R_i, X \\ (D ∪ ⋃_i R_i), D }

returned as a vector `(polytope, obs_label)`.

Labels:
- region `R_i` gets label `i`
- neutral region gets `neutral_obs`
- terminal set `D` gets `terminal_obs`
"""
function build_observation_partition(
    X,
    D,
    regions;
    neutral_obs::Int = 0,
    terminal_obs::Int = -1,
    atol::Float64 = 0.0,
)
    Xh = _as_hpolytope(X)
    Dh = _as_hpolytope(D)
    Rh = [_as_hpolytope(R) for R in regions]

    out = Vector{Tuple{SemiLinearSet, Int}}()

    for (i, R) in enumerate(Rh)
        I = set_intersection(Xh, R)
        if is_nonempty_set(I)
            push!(out, (SemiLinearSet(I), i))
        end
    end

    excluded = SemiLinearSet(vcat([Dh], Rh))
    neutral_set = set_difference_decompose(SemiLinearSet(Xh), excluded; atol = atol)
    if is_nonempty_set(neutral_set)
        push!(out, (neutral_set, neutral_obs))
    end

    Dcap = set_intersection(Xh, Dh)
    if is_nonempty_set(Dcap)
        push!(out, (SemiLinearSet(Dcap), terminal_obs))
    end

    return out
end

# ============================================================
# Main PCLF bisimulation algorithm
# ============================================================

"""
    bisimulation_pclf(f, pclf, Γ, obs_partition; verbose = true, atol = 0.0)

Construct the bisimulation quotient on the lifted product system.
Graph edges are stored as `(source_node, destination_node, mode_label)`.
"""
function bisimulation_pclf(
    f::HybridSystems.HybridSystem,
    pclf::PCLF.PCLF,
    Γ::AbstractVector{<:Real},
    obs_partition::AbstractVector{<:Tuple{<:SemiLinearSet, Int}};
    verbose::Bool = true,
    atol::Float64 = 0.0,
    max_slices::Union{Nothing, Int} = nothing,
)
    A = extract_mode_matrices(f)

    sublevels = build_sublevel_sequence(pclf, Γ)
    slices = build_slice_sequence(sublevels; atol = atol)

    U = typeof(first(pclf.graph.verts))
    SL = SemiLinearSet
    T = PCBisimulationQuotient{SL, U}(slices, obs_partition)

    initialize_partitions!(T)
    initialize_terminal_transitions!(T, pclf)

    incoming_by_dm = group_edges_by_dest_mode(pclf.graph.edges, U)

    N = length(Γ)
    Niter = isnothing(max_slices) ? N : min(max_slices, N)

    refine_count = 0

    for i in 1:Niter
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

"""
    build_sublevel_sequence(pclf, Γ)

Return a dictionary mapping each graph node `s` to the list
`[P^{(s)}_{Γ_1}, ..., P^{(s)}_{Γ_N}]`.
"""
function build_sublevel_sequence(pclf::PCLF.PCLF, Γ::AbstractVector{<:Real})
    U = typeof(first(pclf.graph.verts))
    sublevels = Dict{U, Vector{Poly}}()
    for s in pclf.graph.verts
        piece = pclf.pieces[s]
        sublevels[s] = [_as_hpolytope(PCLF.get_sublevel_set(piece, Float64(γ))) for γ in Γ]
    end
    return sublevels
end

"""
    build_slice_sequence(sublevels; atol = 0.0)

For each node `s`, build
    S_1^{(s)} = P_{Γ_1}^{(s)}
    S_i^{(s)} = P_{Γ_i}^{(s)} \\ P_{Γ_{i-1}}^{(s)},  i>=2

Each slice is stored as a vector of H-polytopes.
"""
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
# Initialization
# ============================================================

"""
    initialize_partitions!(T)

Build the initial node-dependent partitions:
    P_0^(s) = { O ∩ S_i^(s) : O ∈ P_X, S_i^(s) slice, intersection nonempty }.
"""
function initialize_partitions!(T::PCBisimulationQuotient{SemiLinearSet, U}) where {U}
    for (s, slice_list) in T.slices
        for (i, Sset) in enumerate(slice_list)
            for (ObsSet, obs) in T.obs_partition
                I = set_intersection(Sset, ObsSet)
                if is_nonempty_set(I)
                    add_state!(T, s, I, obs, i)
                end
            end
        end
    end
    return T
end

"""
    initialize_terminal_transitions!(T, pclf)

For all terminal states in slice 1, add graph-induced transitions
along every edge `(s,d,m)`.
"""
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

# this could be remove, and adapt the...
@recipe function f(obs_partition::Vector{Tuple{HPolytope, Int}})
    palette = [:gray, :red, :green, :blue, :orange, :purple]

    seen = Set{Int}()

    for (P, obs) in obs_partition
        c = palette[mod1(obs + 2, length(palette))]

        label_str = obs ∈ seen ? "" : "O $obs"
        push!(seen, obs)

        @series begin
            fillalpha := 0.35
            linecolor := c
            fillcolor := c
            label := label_str
            P
        end
    end
end

end # module
