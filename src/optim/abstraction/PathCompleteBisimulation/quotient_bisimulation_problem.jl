
module PathCompleteBisimulation

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

include("geometry_interface.jl")
include("bisimulation_quotient.jl")
include("sublevel_support.jl")

mutable struct OptimizerQuotientBisimulation{T} <: MOI.AbstractOptimizer
    # --- user inputs ---
    quotient_bisimulation_problem::Union{Nothing, PR.QuotientBisimulationProblem}
    pclf::Union{Nothing, PCLF.PCLF}
    Γ::Union{Nothing, Vector{Float64}}
    obs_partition::Union{Nothing, Vector{Tuple{Poly, Int}}}
    verbose::Bool
    atol::T
    num_levels::Union{Nothing, Int}
    tau::Union{Nothing, Float64}

    # --- results ---
    raw_bisimulation::Any
    abstraction_construction_time_sec::T

    function OptimizerQuotientBisimulation{T}() where {T}
        return new{T}(
            nothing,    # quotient_bisimulation_problem
            nothing,    # pclf
            nothing,    # Γ
            nothing,    # obs_partition
            true,       # verbose
            zero(T),    # atol
            nothing,    # num_levels
            nothing,    # tau
            nothing,    # raw_bisimulation
            zero(T),    # solve time
        )
    end
end

OptimizerQuotientBisimulation() = OptimizerQuotientBisimulation{Float64}()

MOI.is_empty(opt::OptimizerQuotientBisimulation) =
    opt.quotient_bisimulation_problem === nothing

function MOI.set(
    model::OptimizerQuotientBisimulation,
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

function MOI.get(model::OptimizerQuotientBisimulation, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec
end

function MOI.get(model::OptimizerQuotientBisimulation, param::MOI.RawOptimizerAttribute)
    name = Symbol(param.name)
    if !hasproperty(model, name)
        error("Unknown optimizer attribute: $(param.name)")
    end
    return getproperty(model, name)
end

function reset!(model::OptimizerQuotientBisimulation)
    model.raw_bisimulation = nothing
    model.slices = nothing
    model.abstraction_construction_time_sec = 0.0
    return model
end

function _validate_model(
    model::OptimizerQuotientBisimulation,
    required_fields::Vector{Symbol},
)
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set `$(field)`. Missing required field in OptimizerQuotientBisimulation.",
            )
        end
    end
end

function MOI.optimize!(opt::OptimizerQuotientBisimulation)
    t_ref = time()

    _validate_model(opt, [:quotient_bisimulation_problem, :pclf])

    prob = opt.quotient_bisimulation_problem

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
    )

    opt.raw_bisimulation = T
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

    out = Tuple{Poly, Int}[]

    for (i, R) in enumerate(Rh)
        I = set_intersection(Xh, R)
        if is_nonempty_set(I)
            push!(out, (I, i))
        end
    end

    excluded = Poly[Dh]
    append!(excluded, Rh)
    neutral_parts = set_difference_decompose(Xh, excluded; atol = atol)
    for P in neutral_parts
        if is_nonempty_set(P)
            push!(out, (P, neutral_obs))
        end
    end

    Dcap = set_intersection(Xh, Dh)
    if is_nonempty_set(Dcap)
        push!(out, (Dcap, terminal_obs))
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
    obs_partition::AbstractVector{<:Tuple{<:Poly, Int}};
    verbose::Bool = true,
    atol::Float64 = 0.0,
)
    A = extract_mode_matrices(f)

    sublevels = build_sublevel_sequence(pclf, Γ)
    slices = build_slice_sequence(sublevels; atol = atol)

    U = typeof(first(pclf.graph.verts))
    T = PCBisimulationQuotient{Poly, U}(slices, obs_partition)

    initialize_partitions!(T)
    initialize_terminal_transitions!(T, pclf)

    N = length(Γ)

    refine_count = 0

    for i in 1:2 # N
        verbose && println("Current slice = $i")

        # Stored order in LabDigraph is (source, destination, label)
        for (s, d, m) in pclf.graph.edges
            verbose && println("  edge = ($s, $m, $d)")

            # Target cells in slice i of destination node d
            target_ids = [
                qid for qid in get(T.part_ids, d, Int[]) if
                haskey(T.states, qid) && T.states[qid].slice == i
            ]
            for qid in target_ids
                haskey(T.states, qid) || continue
                q = T.states[qid]
                preP = preimage_linear(q.set, A[Int(m)])

                # Only refine source cells in strictly outer slices
                source_ids = [
                    pid for pid in get(T.part_ids, s, Int[]) if
                    haskey(T.states, pid) && T.states[pid].slice > i
                ]

                for pid in copy(source_ids)
                    haskey(T.states, pid) || continue
                    # deterministic speed-up: if this state already has a successor for mode m,
                    # do not refine it again for the same mode
                    if any(tr[1] == Int(m) for tr in T.states[pid].next)
                        continue
                    end
                    refine_one_state!(T, pid, preP, Int(m), qid; atol = atol)
                    refine_count += 1

                    if verbose && refine_count % 1000 == 0
                        @info "Refinement progress" refine_count slice=i edge=(s, m, d)
                    end
                end
            end
        end
    end

    verbose && @info "Total number of refinements" refine_count

    return T
end

# ============================================================
# Switched system helpers
# ============================================================

"""
    extract_mode_matrices(f)

Extract the linear reset / mode matrices from a switched `HybridSystem`.
"""
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
    slices = Dict{U, Vector{Vector{Poly}}}()

    for (s, Ps) in sublevels
        Ns = length(Ps)
        local_slices = Vector{Vector{Poly}}(undef, Ns)

        local_slices[1] = [Ps[1]]
        for i in 2:Ns
            local_slices[i] = set_difference_decompose(Ps[i], Ps[i - 1]; atol = atol)
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
function initialize_partitions!(T::PCBisimulationQuotient{Poly, U}) where {Poly, U}
    for (s, slice_list) in T.slices
        for (i, slice_parts) in enumerate(slice_list)
            for Sset in slice_parts
                for (ObsSet, obs) in T.obs_partition
                    I = set_intersection(Sset, ObsSet)
                    if is_nonempty_set(I)
                        add_state!(T, s, I, obs, i)
                    end
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

"""
    refine_one_state!(T, qid, preP, mode, target_qid; atol = 0.0)

Split state `qid` by `preP`.
- difference pieces inherit transitions,
- intersection piece inherits transitions and gets `(mode,target_qid)`.
"""
function refine_one_state!(
    T::PCBisimulationQuotient{Poly, U},
    qid::Int,
    preP::Poly,
    mode::Int,
    target_qid::Int;
    atol::Float64 = 0.0,
) where {U}
    haskey(T.states, qid) || return false
    q = T.states[qid]

    I = set_intersection(q.set, preP)
    if !is_significant_set(I; min_width = atol > 0 ? atol : 1e-8)
        return false
    end

    Dparts = set_difference_decompose(q.set, preP; atol = atol)

    old_next = copy(q.next)
    old_obs = q.obs
    old_node = q.node
    old_slice = q.slice

    remove_state!(T, qid)

    # Difference pieces
    for D in Dparts
        if is_significant_set(D; min_width = atol > 0 ? atol : 1e-8)
            new_id = add_state!(T, old_node, D, old_obs, old_slice)
            T.states[new_id].next = copy(old_next)
        end
    end

    # Intersection piece
    inter_id = add_state!(T, old_node, I, old_obs, old_slice)
    T.states[inter_id].next = copy(old_next)
    add_transition!(T, inter_id, mode, target_qid)

    return true
end

function is_significant_set(P::Poly; min_width::Float64 = 1e-8)
    is_nonempty_set(P) || return false

    # keep only full-dimensional sets if possible
    try
        dim(P) == 2 || return false
    catch
    end

    # crude width check from support in coordinate directions

    try
        xsup = σ([1.0, 0.0], P)
        xinf = -σ([-1.0, 0.0], P)
        ysup = σ([0.0, 1.0], P)
        yinf = -σ([0.0, -1.0], P)

        (xsup - xinf) > min_width || return false
        (ysup - yinf) > min_width || return false
    catch
    end

    return true
end

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
