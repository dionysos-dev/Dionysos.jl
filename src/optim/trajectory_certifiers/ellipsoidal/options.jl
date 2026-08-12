"""
    AdaptiveLinearizationBoxOptions(; kwargs...)

Options of the adaptive linearization-box search (backward direction only): grow the
box on `:lmi_infeasible`, grow to the required radii on `:inconsistent_box`, then
select among box scales according to `objective`.

- `enabled` — turn the adaptive search on (`false` = fixed-box mode using
  `ChainOptions.linearization_δx/δu`);
- `ΔX_initial`/`ΔX_min`/`ΔX_max`, `ΔU_initial`/`ΔU_min`/`ΔU_max` — per-axis
  state/input box radii: the starting box and its clamp bounds (required);
- `growth` — multiplicative box growth on `:lmi_infeasible`;
- `safety` — inflation factor over the required radii on `:inconsistent_box`;
- `max_iters` — cap on adaptive iterations per step;
- `atol` — box-consistency tolerance;
- `search_scales` — box scales tried around the first consistent box
  (`:max_volume` only);
- `objective` — `:first_consistent` (accept the first consistent box — cheapest)
  or `:max_volume` (re-solve on every `search_scales` entry and keep the largest
  certified ellipsoid — ~`length(search_scales)`× the SDPs per step).
"""
Base.@kwdef struct AdaptiveLinearizationBoxOptions
    enabled::Bool = true
    ΔX_initial::Vector{Float64}
    ΔX_min::Vector{Float64}
    ΔX_max::Vector{Float64}
    ΔU_initial::Vector{Float64}
    ΔU_min::Vector{Float64}
    ΔU_max::Vector{Float64}
    growth::Float64 = 1.5
    safety::Float64 = 1.05
    max_iters::Int = 30
    atol::Float64 = 1e-8
    search_scales::Vector{Float64} = [1.0]
    objective::Symbol = :first_consistent
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
- `linearization_δx/δu` — fixed linearization-box radii (fixed mode);
- `adaptive_boxes` — [`AdaptiveLinearizationBoxOptions`](@ref) or `nothing`;
- `objective` — size term of the SDP objective: `:maximin` (largest smallest
  semi-axis — collapse-proof, the default), `:logdet` (true volume), or `:trace`;
- `domain_cap` — make the synthesis domain- and box-aware: cap every funnel
  inside the state domain AND inside the current linearization box by
  construction (per-step SOC slabs, `source_cap` of
  `ST.solve_transition_backward`). Without it a size-maximizing objective grows
  funnels past `X` (state-domain gate rejects a posteriori) or past the box
  (state-side inconsistency drives the adaptive search into ever-bigger boxes
  whose Hessian bounds kill the LMI); with it, box scales become a clean
  line-search dial for the largest certifiable funnel.

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

    # Used when adaptive_boxes === nothing or adaptive_boxes.enabled == false
    linearization_δx::Vector{Float64} = Float64[]
    linearization_δu::Vector{Float64} = Float64[]

    adaptive_boxes::Union{Nothing, AdaptiveLinearizationBoxOptions} = nothing
    objective::Symbol = :maximin
    # :vertices enumerates the 2ⁿ Lipschitz-box corners (exact, tightest — measured
    # +33% certified runway on the double pendulum vs :ball); :ball wraps them in
    # one scalar norm ball (constant block count — cheapest at n ≥ 4); :john_ball
    # wraps them in the box's John ellipsoid √n·diag(Lip)·B (per-axis radii at
    # :ball's block count — the scalable middle ground).
    remainder_model::Symbol = :vertices
    domain_cap::Bool = false
    # Retry a failed step as ONE two-step transition into E_{k+2} through the
    # already-synthesized κ_{k+1}, skipping the intermediate containment. Across
    # a rescued pair the funnel guarantee is the composition E_k → (2 steps) →
    # E_{k+2}; E_{k+1} membership is not claimed. Requires adaptive_boxes.
    # MEASURED NEGATIVE on the double pendulum (kept as
    # machinery): by the first failure the downstream funnels are needles whose
    # margin the composition's ℓ-offset and constant e₂ budget already consume,
    # and a stride-2 chain from the terminal is strictly worse than 1-step
    # (fixed mid-κ + worst-case e₂ lose more than the skipped containment
    # saves). The true 2-step coupling advantage needs joint bilinear gain
    # synthesis — nonconvex, out of scope.
    two_step_rescue::Bool = false

    r_min::Float64 = 0.0
    check_state_domain::Bool = true
    check_terminal::Bool = true
end

"""
    ForwardOptions(; kwargs...)

Options of the forward certification chain.

- `target_mode` — `:free` (target shape `Q₂` is a decision variable, trace
  objective, conditioning sandwich `q_min`/`q_max`) or `:fixed` (shape follows the
  entry ellipsoid's, only the scale `α` is free — its per-step value is the
  contraction profile);
- `α_max` — fail-fast guard: a step whose tube scale — `tr(Q_{k+1})/tr(Q₁)`,
  relative to the ENTRY shape in both target modes — exceeds this is rejected
  (`Inf` disables);
- `entry_shape` — LazySets shape matrix of the entry ellipsoid `E₁`, or `nothing`
  to circumscribe the problem's initial set (centered at the trajectory start);
- `maxδu`, `λ`, `transition_cost`, `linearization_δu` as in [`ChainOptions`](@ref).
  The two target modes want OPPOSITE λ regimes (measured on a linear system):
  in `:fixed` mode α is dimensionless, and λ ≪ 1 polishes the min-α solution
  onto the strict-PSD boundary where the solver tolerance eats the ε margin and
  the a-posteriori validation rejects a "solved" step — use λ ≈ 0.5; in `:free`
  mode the trace term is in absolute tube units, and a large λ lets the cost
  term buy input effort by drifting the tube off the nominal until input
  feasibility dies — keep λ small (the 0.01 default) and instead set `q_min`
  near the entry scale, since at the loose default floor the trace objective
  needle-collapses the tube and the next step dies of source conditioning;
- `linearization_δx_margin ≥ 1` inflates the (known) state box handed to the
  Hessian bound, buying u-side slack;
- `remainder_model` — `:vertices` or `:ball` (`:john_ball` is backward-only);
- gates: `r_min`, `check_state_domain`, `check_terminal` as in [`ChainOptions`](@ref).
"""
Base.@kwdef mutable struct ForwardOptions
    target_mode::Symbol = :free
    α_max::Float64 = Inf
    q_min::Float64 = 1e-9
    q_max::Float64 = 1e9
    entry_shape::Union{Nothing, Matrix{Float64}} = nothing
    maxδu::Float64 = 0.5
    λ::Float64 = 0.01
    transition_cost::Union{Nothing, UT.QuadraticStateControlFunction, Matrix{Float64}} =
        nothing
    linearization_δu::Vector{Float64} = Float64[]
    linearization_δx_margin::Float64 = 1.1
    r_min::Float64 = 0.0
    check_state_domain::Bool = true
    check_terminal::Bool = true
    remainder_model::Symbol = :vertices
end
