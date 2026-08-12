# Forward certification (plan.md §4.4): the tube is *propagated* from a given entry
# ellipsoid, one `ST.solve_transition_forward` SDP per step. The source of every
# step is data, so the linearization box is known exactly before solving — the
# x-side consistency gate holds by construction and no adaptive boxes exist here.
# Failure is late but informative (the contraction profile shows where the tube
# inflates); the terminal gate `E_{K+1} ⊆ target_set` closes the specification.
#
# `ForwardOptions` lives in options.jl; the certifier front in certifier.jl.

# Circumscribed ellipsoid of the initial set, centered at the trajectory start:
# for a box with vertices vⱼ, the diagonal shape with semi-axes √n·maxⱼ|vⱼᵢ − cᵢ|
# contains every vertex (Cauchy–Schwarz), hence the box.
function _entry_ellipsoid(problem, x1, opts::ForwardOptions)
    if opts.entry_shape !== nothing
        return LazySets.Ellipsoid(collect(float.(x1)), Matrix{Float64}(opts.entry_shape)),
        nothing
    end
    verts = try
        LazySets.vertices_list(problem.initial_set)
    catch
        return nothing,
        "cannot circumscribe the initial set (no vertex list) — set entry_shape"
    end
    isempty(verts) && return nothing, "empty initial set"
    n = length(x1)
    h = [maximum(abs(v[i] - x1[i]) for v in verts) for i in 1:n]
    all(h .>= 0) || return nothing, "degenerate initial set"
    a = sqrt(n) .* max.(h, 1e-9)
    return LazySets.Ellipsoid(collect(float.(x1)), Matrix(LA.Diagonal(a .^ 2))), nothing
end

function forward_step!(ctx::ChainContext, k::Int, E_k, Qhat, Qref)
    opts = ctx.options
    @assert !isempty(opts.linearization_δu) "Set options.linearization_δu."

    xk = collect(LazySets.center(E_k))
    uk = collect(ctx.us[k])
    xnext = ctx.xs[k + 1]
    # The state box is known exactly from E_k — the whole point of forward mode.
    δx = opts.linearization_δx_margin .* _ellipsoid_axis_radii(E_k)
    δu = collect(Float64, opts.linearization_δu)

    approx = ST.build_affine_approximation(ctx.affine_provider, xk, uk; δx = δx, δu = δu)

    result = ST.solve_transition_forward(
        approx.system,
        E_k,
        collect(float.(xnext)),
        uk,
        approx.Uformat,
        approx.Wformat,
        ctx.S,
        approx.lipschitz,
        ctx.backend;
        target_shape = opts.target_mode === :fixed ? Qhat : nothing,
        maxδu = opts.maxδu,
        λ = opts.λ,
        q_min = opts.q_min,
        q_max = opts.q_max,
        remainder_model = opts.remainder_model,
    )

    if !result.feasible
        return StepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (; Xbar_radius = δx, Ubar_radius = δu, step_status = :forward_infeasible),
        )
    end

    # u-side consistency (the x-side holds by construction): the controller image
    # over E_k must stay inside the box the Hessian bound was taken on.
    uc, ru = _controller_image_axis_radii(result.controller, E_k)
    required_δu = abs.(uc .- collect(Float64, uk)) .+ ru
    if !all(required_δu .<= δu .+ 1e-8)
        return StepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (;
                Xbar_radius = δx,
                Ubar_radius = δu,
                required_U_radius = required_δu,
                step_status = :forward_inconsistent_δu,
            ),
        )
    end

    E_next = result.target
    scale = LA.tr(Matrix{Float64}(LazySets.shape_matrix(E_next))) / LA.tr(Qref)

    return StepRecord(
        k,
        :ok,
        result.cost,
        E_next,
        result.controller,
        (;
            Xbar_radius = δx,
            Ubar_radius = δu,
            required_U_radius = required_δu,
            contraction = scale,
            L = approx.lipschitz,
        ),
    )
end

function _run_forward_chain(ctx::ChainContext)
    t0 = time()
    problem = ctx.problem
    opts = ctx.options

    E_1, entry_reason = _entry_ellipsoid(problem, ctx.xs[1], opts)
    if E_1 !== nothing && entry_reason === nothing
        # The entry enters no StepRecord — gate it here or never.
        entry_reason = endpoint_gate(
            E_1,
            problem.system.X,
            opts.r_min,
            opts.check_state_domain,
            "entry",
        )
        entry_reason === nothing || (E_1 = nothing)
    end
    E_1 === nothing && return _failed_before_start(0, t0, nothing, entry_reason)

    Qref = Matrix{Float64}(LazySets.shape_matrix(E_1))
    Qhat = opts.target_mode === :fixed ? Qref : nothing

    # Real coverage, not an assumption: with a user-supplied `entry_shape` the
    # entry ellipsoid has no relation to the initial set.
    initial_cov = initial_coverage(problem.initial_set, E_1)

    domain_checked = opts.check_state_domain && _state_domain_supported(problem.system.X)

    steps = StepRecord[]
    ellipsoids = [E_1]
    E_k = E_1

    for k in 1:(ctx.K)
        rec = forward_step!(ctx, k, E_k, Qhat, Qref)
        rec = apply_gates(rec, problem, opts)
        push!(steps, rec)

        if rec.status != :ok
            return _assemble_result(
                false,
                k,
                t0,
                steps,
                copy(ellipsoids),
                nothing,
                initial_cov,
                domain_checked,
            )
        end

        E_k = rec.ellipsoid
        push!(ellipsoids, E_k)
    end

    terminal_contained =
        opts.check_terminal ? terminal_containment(E_k, problem.target_set) : nothing
    success = terminal_contained !== false

    return _assemble_result(
        success,
        success ? nothing : ctx.K + 1,
        t0,
        steps,
        copy(ellipsoids),
        terminal_contained,
        initial_cov,
        domain_checked,
    )
end
