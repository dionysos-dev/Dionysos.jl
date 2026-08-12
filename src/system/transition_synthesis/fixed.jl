# Fixed-shape kernel: both ellipsoids given, decision variables (K, ℓ).

# The target enters the reachability blocks as its inverse quadratic form,
# which is exactly the LazySets shape matrix Q2 — no inversion needed.
function _fixed_shape_kernel(A, B, c, Wcols, U, Λ, c1, P1, c2, Q2, sdp_solver)
    n = length(c)
    Nw = size(Wcols, 2)
    Nu = length(U)
    m = size(U[1], 2)

    km = _KernelModel(sdp_solver, 1e-10)
    model = km.model
    @variable(model, K[i = 1:m, j = 1:n])
    @variable(model, ell[i = 1:m, j = 1:1])
    @variable(model, beta[i = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)
    _require_nonneg!(km, beta, tau)

    @expressions(model, begin
        At, A + B * K
        ct, c + B * ell
    end)

    z = zeros(n, 1)

    for i in 1:Nw
        aux = @expression(model, A * hcat(c1) + hcat(ct) - hcat(c2) + hcat(Wcols[:, i]))
        _pose_psd!(km, _reach_block, beta[i], P1, At, aux, Q2, n)
    end

    _pose_inputs!(km, tau, U, P1, K, ell, n)

    G = [LA.I transpose(K) z]
    d = [transpose(c1) transpose(ell) 1]
    _pose_cost!(km, γ, J, P1, G, d, Λ, n)

    @objective(model, Min, J)

    optimize!(model)

    if !_solver_accepted(model)
        @debug "solve_transition: infeasible SDP" termination_status(model)
        return false, nothing, nothing
    end
    _validated(km, "solve_transition") || return false, nothing, nothing

    return true, [_solved(K) _solved(ell)], value(J)
end

"""
    solve_transition(affsys, source, target, U, W, cost, sdp_solver) -> TransitionResult

Synthesize an affine controller `u(x) = K(x − c₁) + ℓ` certifying that every
state of the `source` ellipsoid reaches the `target` ellipsoid in one step of
the affine system `affsys` (`x⁺ = Ax + Bu + c (+ Dw)`), for every disturbance
vertex, while minimizing an upper bound on the worst-case transition cost
`‖cost · [x; u; 1]‖²`.

# Arguments
- `affsys`: an affine system (`(Noisy)ConstrainedAffineControlDiscreteSystem`
  or `ConstrainedAffineControlMap`);
- `source`, `target`: `LazySets.Ellipsoid`s;
- `U`: input constraints — a set (`LazySets.AbstractHyperrectangle`, `LazySets.Ellipsoid`,
  `LazySets.IntersectionArray`) or a preformatted list (`format_input_set`);
- `W`: disturbance vertices — one per column (at least one; a single zero vertex
  for deterministic dynamics), or a box (`LazySets.AbstractHyperrectangle`);
  mapped through the system's noise matrix `D` when the system has one;
- `cost`: the PSD cost matrix `S` of size `(nx+nu+1) × (nx+nu+1)` bounding
  `[x; u; 1]ᵀ·S·[x; u; 1]`, or a `UT.QuadraticStateControlFunction`; the factor
  the LMIs consume is taken internally (`_cost_factor`);
- `sdp_solver`: a JuMP-compatible SDP optimizer (Clarabel, Mosek, …).

Every returned certificate is re-validated numerically at the solver's solution —
a transition reported feasible is PSD-certified independently of the solver's
termination status.
"""
function solve_transition(
    affsys::AffineSys,
    source::LazySets.Ellipsoid,
    target::LazySets.Ellipsoid,
    U,
    W,
    cost,
    sdp_solver,
)
    c1 = LazySets.center(source)
    feasible, kappa, J = _fixed_shape_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _check_noise_vertices(_state_noise(affsys, W)),
        _input_matrices(U),
        _cost_factor(cost),
        c1,
        UT.get_quadratic_form(source),
        LazySets.center(target),
        LazySets.shape_matrix(target),
        sdp_solver,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    return TransitionResult(true, MS.AffineMap(K, ℓ - K * c1), J, nothing)
end
