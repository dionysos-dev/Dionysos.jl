# The chain context: everything `run_chain!` needs, built once per certification —
# problem, trajectory data, affine-approximation provider, transition cost, SDP
# backend, options — plus the backward transition-solve wrapper.
#
# Per-step state scaling used to live here; it was superseded by globally
# normalized dynamics (`ST.normalized_symbolic_provider`), which keeps the
# conditioning benefit with EXACT Hessian bounds in the working frame instead of
# the scaling's ~1/σ_min(T)² remainder tax.

struct ChainContext{P, AP, TX, TU, TS, TB, OPTS}
    problem::P
    affine_provider::AP
    xs::TX
    us::TU
    K::Int
    S::TS
    backend::TB
    options::OPTS
end

function _identity_transition_cost(nx::Int, nu::Int)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

function _transition_cost_matrix(cost, nx, nu)
    S =
        cost === nothing ? _identity_transition_cost(nx, nu) :
        Matrix{Float64}(UT._cost_matrix(cost))
    @assert size(S) == (nx + nu + 1, nx + nu + 1)
    return S
end

# `options` is ChainOptions (backward) or ForwardOptions (forward) — the context
# only reads `transition_cost` here; the chains read their own fields.
function build_context(
    problem::PR.ProblemType,
    traj::ST.Trajectory,
    affine_provider,
    backend,
    options,
)
    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))

    K = length(us)
    @assert length(xs) == K + 1 "Expected length(xs) == length(us) + 1."

    nx = length(xs[1])
    nu = length(us[1])
    S = _transition_cost_matrix(options.transition_cost, nx, nu)

    return ChainContext(problem, affine_provider, xs, us, K, S, backend, options)
end

# ------------------------------------------------------------
# Transition solve
# ------------------------------------------------------------

# Per-step slab cap from the state domain: the distance from the nominal center to
# each domain face, shaved so solver-tolerance solutions still pass the state-domain
# gate. `nothing` when the domain is unsupported or the nominal sits outside it (the
# gate then reports the violation as before).
function _domain_cap(X, xk)
    box = X isa UT.SetMinus ? UT.minus_included(X) : X
    box isa LazySets.AbstractHyperrectangle || return nothing
    cap = min.(LazySets.high(box) .- xk, xk .- LazySets.low(box))
    all(cap .> 0) || return nothing
    return cap .* (1 - 1e-6)
end

# Merge the domain cap with the linearization-box cap into the kernel's `source_cap`
# slab rows; the shave keeps solver-tolerance solutions inside the gates.
function _merged_source_cap(ctx::ChainContext, xk, box_cap)
    source_cap =
        ctx.options.domain_cap ? _domain_cap(ctx.problem.system.X, collect(Float64, xk)) :
        nothing
    if box_cap !== nothing
        shaved = collect(Float64, box_cap) .* (1 - 1e-6)
        source_cap = source_cap === nothing ? shaved : min.(source_cap, shaved)
    end
    return source_cap
end

function _solve_transition(ctx::ChainContext, approx, E_next, xk, uk; box_cap = nothing)
    # Cap the source inside the state domain AND the linearization box: state-side
    # box consistency then holds by construction, so the adaptive search
    # line-searches box scales for the largest certifiable funnel instead of
    # chasing a size-maximizing SDP with ever-bigger boxes (whose Hessian bounds
    # eventually kill the LMI — measured `lmi_infeasible_at_max_box` on the double
    # pendulum's mid-swing).
    result = ST.solve_transition_backward(
        approx.system,
        E_next,
        xk,
        uk,
        approx.Uformat,
        approx.Wformat,
        ctx.S,
        approx.lipschitz,
        ctx.backend;
        maxδx = ctx.options.maxδx,
        maxδu = ctx.options.maxδu,
        λ = ctx.options.λ,
        objective = ctx.options.objective,
        remainder_model = ctx.options.remainder_model,
        source_cap = _merged_source_cap(ctx, xk, box_cap),
    )
    return result.source, result.controller, result.cost
end
