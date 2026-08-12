# Free-target kernel (forward direction): the source ellipsoid is data — so the
# linearization deviation δx is a *constant* and the remainder is known before
# solving — and the target enters the reach blocks only as its shape matrix Q₂,
# where it is linear. A free target is the naturally convex direction of ellipsoid
# calculus: no L-congruence, no volume cone, no adaptive boxes (plan.md §4.4).

function _free_target_kernel(
    A,
    B,
    c,
    Wcols,
    U,
    Λ,
    c1,
    P1,
    δx_const,
    u_ref,
    c2,
    Qhat,          # target shape matrix for the α-mode; `nothing` = free-shape mode
    Lip,
    sdp_solver;
    maxδu,
    λ,
    q_min,
    q_max,
    remainder_model::Symbol,
)
    remainder_model === :john_ball &&
        error("remainder_model = :john_ball is implemented for the backward kernel only.")
    nx = length(c)
    nu = size(U[1], 2)
    spec = _remainder_spec(remainder_model, Lip, nx)
    νs = [Wcols[:, i] for i in 1:size(Wcols, 2)]
    Nx = spec.ball_like ? 1 : length(spec.vertices)
    Nw = length(νs)
    Nu = length(U)

    km = _KernelModel(sdp_solver, 1e-8)
    model = km.model
    @variable(model, K[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx, j = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)
    muball = spec.ball_like ? @variable(model, [j = 1:Nw], lower_bound = 0.0) : nothing
    _require_nonneg!(km, beta, tau)
    muball === nothing || _require_nonneg!(km, muball)

    local α, Q2
    if Qhat !== nothing
        @variable(model, α >= 0)
        Q2 = α .* Qhat
    else
        @variable(model, Q2[i = 1:nx, j = 1:nx], PSD)
        # Conditioning sandwich: thin targets are fine for the certificate but the
        # next step inverts Q₂ into its source form.
        @constraint(model, Q2 >= q_min * _eye(nx), PSDCone())
        @constraint(model, q_max * _eye(nx) >= Q2, PSDCone())
    end

    @expressions(model, begin
        At, A + B * K
        ct, c + B * ell
    end)

    z = zeros(nx, 1)
    base_aux = @expression(model, A * hcat(c1) + hcat(ct) - hcat(c2))
    scale = @expression(model, δx_const + δu)
    _pose_reach!(km, spec, beta, muball, P1, At, base_aux, νs, Q2, scale, nx)

    _pose_inputs!(km, tau, U, P1, K, ell, nx)

    G = [LA.I transpose(K) z]
    d = [transpose(c1) transpose(ell) 1]
    _pose_cost!(km, γ, J, P1, G, d, Λ, nx)

    _pose_proximity!(km, ϕ, P1, K, ell, u_ref, δu, nx, nu)

    @constraint(model, δu <= maxδu^2)
    _check!(() -> _solved(δu) <= maxδu^2 + _VALIDATION_TOL, km)

    size_term = Qhat !== nothing ? α : sum(Q2[i, i] for i in 1:nx)
    @objective(model, Min, λ * J + (1.0 - λ) * size_term)

    optimize!(model)

    if !_solver_accepted(model)
        @debug "solve_transition_forward: infeasible SDP" termination_status(model)
        return false, nothing, nothing, nothing
    end
    _validated(km, "solve_transition_forward") || return false, nothing, nothing, nothing

    Q2v = Qhat !== nothing ? value(α) .* Qhat : _solved(Q2)
    return true, Q2v, [_solved(K) _solved(ell)], value(J)
end

"""
    solve_transition_forward(affsys, source, target_center, u_ref, U, W, cost,
                             lipschitz, sdp_solver;
                             target_shape = nothing, maxδu = 20.0, λ = 0.01,
                             q_min = 1e-9, q_max = 1e9,
                             remainder_model = :vertices) -> TransitionResult

Synthesize an affine controller *and the smallest target ellipsoid* centered at
`target_center` that the controller certifiably reaches from the **given** `source`
ellipsoid in one step of `affsys`, robustly to the disturbance vertices and to the
linearization error bounded by `lipschitz` — whose state deviation is the *known*
`λ_max(Q₁)` of the source (a SQUARED semi-axis: the same squared-radius
convention as the backward kernels' `δx`), so the remainder is fixed data (no
adaptive boxes).

Target modes: with `target_shape` (a LazySets shape matrix `Q̂`) only the scale `α`
is free (`min λ·J + (1−λ)·α` — `α` is the per-step contraction number); without it
the shape `Q₂` itself is a decision variable (`min λ·J + (1−λ)·tr(Q₂)`) inside the
conditioning sandwich `q_min·I ⪯ Q₂ ⪯ q_max·I`. Both are single convex SDPs — the
target enters the reach blocks linearly. `remainder_model` accepts `:vertices`
or `:ball` (`:john_ball` is backward-only and errors here). Certificates are
validated at the returned solution as in [`solve_transition`](@ref). Returns a
[`TransitionResult`](@ref) whose `target` is the synthesized ellipsoid.
"""
function solve_transition_forward(
    affsys::AffineSys,
    source::LazySets.Ellipsoid,
    target_center,
    u_ref,
    U,
    W,
    cost,
    lipschitz,
    sdp_solver;
    target_shape = nothing,
    maxδu = 20.0,
    λ = 0.01,
    q_min = 1e-9,
    q_max = 1e9,
    remainder_model::Symbol = :vertices,
)
    c1 = LazySets.center(source)
    P1 = Matrix{Float64}(UT.get_quadratic_form(source))
    δx_const = LA.eigmax(LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(source))))

    feasible, Q2, kappa, J = _free_target_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _check_noise_vertices(_state_noise(affsys, W)),
        _input_matrices(U),
        _cost_factor(cost),
        c1,
        P1,
        δx_const,
        u_ref,
        target_center,
        target_shape,
        lipschitz,
        sdp_solver;
        maxδu = maxδu,
        λ = λ,
        q_min = q_min,
        q_max = q_max,
        remainder_model = remainder_model,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    controller = MS.AffineMap(K, ℓ - K * c1)
    target = LazySets.Ellipsoid(
        collect(float.(target_center)),
        UT._symmetrize(Matrix{Float64}(Q2));
        check_posdef = false,
    )
    return TransitionResult(true, controller, J, nothing, target)
end
