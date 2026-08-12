# Free-source kernel: target given, the source shape (via its square root L)
# is a decision variable; the source center c₁ and a reference input are data.

function _free_source_kernel(
    A,
    B,
    c,
    Wcols,
    U,
    Λ,
    c1,
    u_ref,
    c2,
    Q2,
    Lip,
    sdp_solver;
    maxδx,
    maxδu,
    λ,
    objective::Symbol,
    remainder_model::Symbol,
    source_cap = nothing,
)
    nx = length(c)
    nu = size(U[1], 2)
    spec = _remainder_spec(remainder_model, Lip, nx)
    νs = [Wcols[:, i] for i in 1:size(Wcols, 2)]
    Nx = spec.ball_like ? 1 : length(spec.vertices)
    Nw = length(νs)
    Nu = length(U)

    km = _KernelModel(sdp_solver, 1e-8)
    model = km.model
    @variable(model, L[i = 1:nx, j = 1:nx], PSD)
    @variable(model, F[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx, j = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, δx >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)
    muball = spec.ball_like ? @variable(model, [j = 1:Nw], lower_bound = 0.0) : nothing
    _require_nonneg!(km, beta, tau)
    muball === nothing || _require_nonneg!(km, muball)

    _pose_source_cap!(km, source_cap, L, nx)

    @expressions(model, begin
        At, A * L + B * F
        ct, c + B * ell
    end)

    z = zeros(nx, 1)
    base_aux = @expression(model, A * hcat(c1) + hcat(ct) - hcat(c2))
    scale = @expression(model, δx + δu)
    _pose_reach!(km, spec, beta, muball, _eye(nx), At, base_aux, νs, Q2, scale, nx)

    _pose_inputs!(km, tau, U, _eye(nx), F, ell, nx)

    G = [transpose(L) transpose(F) z]
    d = [transpose(c1) transpose(ell) 1]
    _pose_cost!(km, γ, J, _eye(nx), G, d, Λ, nx)

    _pose_proximity!(km, ϕ, _eye(nx), F, ell, u_ref, δu, nx, nu)
    _pose_psd!(km, _source_radius_block, L, δx, nx)
    _pose_deviation_caps!(km, δx, maxδx, δu, maxδu)

    _pose_source_size_objective!(km, objective, λ, L, J, nx)

    optimize!(model)

    if !_solver_accepted(model)
        @debug "solve_transition_backward: infeasible SDP" termination_status(model)
        return false, nothing, nothing, nothing
    end
    _validated(km, "solve_transition_backward") || return false, nothing, nothing, nothing

    Q1, kappa = _source_solution(L, F, ell, "solve_transition_backward")
    Q1 === nothing && return false, nothing, nothing, nothing
    return true, Q1, kappa, value(J)
end

"""
    solve_transition_backward(affsys, target, source_center, u_ref, U, W, cost,
                              lipschitz, sdp_solver;
                              maxδx = 100.0, maxδu = 20.0, λ = 0.01,
                              objective = :logdet, remainder_model = :vertices,
                              source_cap = nothing) -> TransitionResult

Synthesize an affine controller *and the largest source ellipsoid* centered at
`source_center` that the controller certifiably drives into the `target`
ellipsoid in one step of `affsys`, robustly to the disturbance vertices `W` and
to the linearization error bounded by `lipschitz` (the Lipschitz radii of
`[x; u; 1]`, scaled by the synthesized deviations `δx + δu` — NOTE that `δx` /
`δu` bound the SQUARED state/input deviations, so `maxδx` / `maxδu` are radii
at this boundary while the remainder is quadratic in them).

The `objective` trades the transition-cost bound against the source size
(`min λ·J − (1−λ)·size`): `:maximin` maximizes the smallest semi-axis
(collapse-proof and cone-free — the ellipsoidal certifier's default), `:logdet`
the true volume (this function's default), `:trace` its stable proxy.
`remainder_model` selects how the linearization error enters the reach blocks:
`:vertices` (exact 2ⁿ corners), `:ball` (one scalar norm ball), or `:john_ball`
(the box's John ellipsoid — per-axis radii at `:ball`'s block count).
`u_ref` is the reference input the controller stays `δu`-close to. `source_cap`
(a per-state vector `d`) additionally confines the source to the axis-aligned
slab `|xᵢ − c₁ᵢ| ≤ dᵢ` — one SOC row per state on the shape factor — so a chain
can keep its funnels inside a state domain by construction. Other arguments as
in [`solve_transition`](@ref). Returns a [`TransitionResult`](@ref) whose
`source` is the synthesized ellipsoid.
"""
function solve_transition_backward(
    affsys::AffineSys,
    target::LazySets.Ellipsoid,
    source_center,
    u_ref,
    U,
    W,
    cost,
    lipschitz,
    sdp_solver;
    maxδx = 100.0,
    maxδu = 20.0,
    λ = 0.01,
    objective::Symbol = :logdet,
    remainder_model::Symbol = :vertices,
    source_cap = nothing,
)
    feasible, Q1, kappa, J = _free_source_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _check_noise_vertices(_state_noise(affsys, W)),
        _input_matrices(U),
        _cost_factor(cost),
        source_center,
        u_ref,
        LazySets.center(target),
        LazySets.shape_matrix(target),
        lipschitz,
        sdp_solver;
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
        objective = objective,
        remainder_model = remainder_model,
        source_cap = source_cap,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    controller = MS.AffineMap(K, ℓ - K * source_center)
    source = LazySets.Ellipsoid(
        collect(float.(source_center)),
        UT._symmetrize(Q1);
        check_posdef = false,
    )
    return TransitionResult(true, controller, J, source)
end
