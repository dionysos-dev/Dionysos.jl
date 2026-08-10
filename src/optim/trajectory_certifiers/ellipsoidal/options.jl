"""
    AdaptiveLinearizationBoxOptions

Options of the adaptive linearization-box search (backward direction only): grow the
box on `:lmi_infeasible`, grow to the required radii on `:inconsistent_box`, then
optionally line-search box scales for the largest certified ellipsoid.
"""
struct AdaptiveLinearizationBoxOptions
    enabled::Bool
    ΔX_initial::Vector{Float64}
    ΔX_min::Vector{Float64}
    ΔX_max::Vector{Float64}
    ΔU_initial::Vector{Float64}
    ΔU_min::Vector{Float64}
    ΔU_max::Vector{Float64}
    growth::Float64
    safety::Float64
    max_iters::Int
    atol::Float64
    verbose::Bool
    search_scales::Vector{Float64}
    objective::Symbol
    keep_first_consistent::Bool
end

"""
    ChainOptions(; kwargs...)

Options of the ellipsoidal certification chain.

- `maxδx`, `maxδu` — caps on the synthesized state/input deviations (passed to
  `ST.solve_transition_backward`);
- `λ` — cost-vs-volume trade-off of the backward SDP objective
  (`min λ·J − (1−λ)·volume`); the default `0.01` is volume-dominated — `λ = 1.0`
  would remove the volume incentive entirely and let funnels collapse;
- `terminal_shape` — LazySets shape matrix `Q` of the terminal ellipsoid, or
  `nothing` to inscribe an ellipsoid in the problem's target set around the
  trajectory endpoint, shrunk by `terminal_shrink`;
- `transition_cost` — PSD `[x; u; 1]` cost matrix or `UT.QuadraticStateControlFunction`
  (identity when `nothing`);
- `state_scaling` — optional per-step change of coordinates (see `_scaled_lipschitz`
  for its soundness caveat; prefer globally normalized dynamics);
- `linearization_δx/δu` — fixed linearization-box radii (fixed mode);
- `adaptive_boxes` — [`AdaptiveLinearizationBoxOptions`](@ref) or `nothing`;
- `objective` — size term of the SDP objective: `:maximin` (largest smallest
  semi-axis — collapse-proof, the default), `:logdet` (true volume), or `:trace`.

Soundness gates (plan.md §4.2):

- `r_min` — minimum admissible semi-axis of every funnel ellipsoid (collapse gate;
  `0` disables);
- `check_state_domain` — require every funnel ellipsoid inside the system domain `X`
  and provably disjoint from its holes (reach-avoid gate);
- `check_terminal` — require the terminal ellipsoid inside the problem's target set.
"""
Base.@kwdef mutable struct ChainOptions
    maxδx::Float64 = 0.2
    maxδu::Float64 = 0.5
    λ::Float64 = 0.01
    terminal_shape::Union{Nothing, Matrix{Float64}} = nothing
    terminal_shrink::Float64 = 0.9
    transition_cost::Union{Nothing, UT.QuadraticStateControlFunction, Matrix{Float64}} =
        nothing
    state_scaling::Union{Nothing, Matrix{Float64}} = nothing

    # Used when adaptive_boxes === nothing or adaptive_boxes.enabled == false
    linearization_δx::Vector{Float64} = Float64[]
    linearization_δu::Vector{Float64} = Float64[]

    adaptive_boxes::Union{Nothing, AdaptiveLinearizationBoxOptions} = nothing
    objective::Symbol = :maximin

    r_min::Float64 = 0.0
    check_state_domain::Bool = true
    check_terminal::Bool = true
end
