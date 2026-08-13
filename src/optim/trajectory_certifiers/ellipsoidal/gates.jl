# Soundness gates (plan.md §4.2). Pure functions returning `nothing` (pass) or a
# `String` reason (fail); the chain converts any failure into an :infeasible record.
# A step that clears its SDP but fails a gate is NOT certified — the gates are what
# connect the per-step LMI certificate to the problem's specification.

# ------------------------------------------------------------
# Ellipsoid-vs-box primitives. For an ellipsoid with shape matrix Q the support in
# the axis direction eᵢ is cᵢ + √Qᵢᵢ, so containment in an axis-aligned box is
# exact; disjointness uses the bounding box (sufficient, hence sound).
# ------------------------------------------------------------

function _ellipsoid_bbox(E)
    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    r = sqrt.(max.(0.0, LA.diag(Q)))
    return LazySets.Hyperrectangle(collect(Float64, LazySets.center(E)), r)
end

# Per-axis support excess over the box faces (> 0 on the violated axes).
function _box_excess(E, H::LazySets.AbstractHyperrectangle)
    c = collect(Float64, LazySets.center(E))
    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    cH = collect(Float64, LazySets.center(H))
    rH = [LazySets.radius_hyperrectangle(H, i) for i in 1:LazySets.dim(H)]
    return abs.(c .- cH) .+ sqrt.(max.(0.0, LA.diag(Q))) .- rH
end

_ellipsoid_in_box(E, H::LazySets.AbstractHyperrectangle) = all(_box_excess(E, H) .<= 1e-12)

# Failure text naming the violated axes — a domain rejection without the axes is
# undebuggable from the outside (the gate drops the uncertified ellipsoid).
function _box_violation_reason(E, H)
    excess = _box_excess(E, H)
    axes = findall(excess .> 1e-12)
    return "is not contained in the state domain (axis " *
           "$(join(axes, ", ")): support excess $(round.(excess[axes]; sigdigits = 3)))"
end

_provably_disjoint(E, S) = UT.is_disjoint(_ellipsoid_bbox(E), S)

_each_member(S::LazySets.UnionSetArray) = LazySets.array(S)
_each_member(S) = (S,)

# ------------------------------------------------------------
# Gate A — box consistency (mandatory). An ellipsoid (or controller image) larger
# than the linearization box invalidates the Hessian remainder bound it was
# certified with; the fixed-box mode used to record this without enforcing it.
# ------------------------------------------------------------

function box_consistency_gate(rec::StepRecord)
    rec.status == :ok || return nothing
    s = rec.summary
    (s.required_X_radius === nothing || s.required_U_radius === nothing) && return nothing
    # A zero Hessian bound (exactly affine dynamics) makes the remainder zero for
    # any box — the box then carries no soundness content and cannot be violated.
    s.L !== nothing && all(abs.(collect(Float64, s.L)) .<= 1e-14) && return nothing
    _box_contains_required_radii(
        s.Xbar_radius,
        s.Ubar_radius,
        s.required_X_radius,
        s.required_U_radius;
        atol = 1e-8,
    ) && return nothing
    return "certified ellipsoid exceeds its linearization box " *
           "(required δx = $(round.(s.required_X_radius; sigdigits = 3)) vs box " *
           "$(round.(s.Xbar_radius; sigdigits = 3)); required δu = " *
           "$(round.(s.required_U_radius; sigdigits = 3)) vs box " *
           "$(round.(s.Ubar_radius; sigdigits = 3)))"
end

# ------------------------------------------------------------
# Gate G — collapse floor: every funnel ellipsoid must keep a minimum semi-axis.
# ------------------------------------------------------------

function collapse_gate(rec::StepRecord, r_min::Float64)
    rec.status == :ok || return nothing
    r_min <= 0.0 && return nothing
    Q = LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(rec.ellipsoid)))
    rmin = sqrt(max(0.0, LA.eigmin(Q)))
    rmin >= r_min && return nothing
    return "funnel collapsed: min semi-axis $(round(rmin; sigdigits = 3)) < r_min $(r_min)"
end

# ------------------------------------------------------------
# Gate C — state domain / reach-avoid: every funnel ellipsoid inside the domain and
# provably disjoint from its holes (obstacles). Supports plain boxes and
# `UT.set_minus` domains; other domain types are not checked (reported once by the
# chain as `state_domain_checked = false`).
# ------------------------------------------------------------

_state_domain_supported(X) = X isa LazySets.AbstractHyperrectangle || X isa UT.SetMinus

# Ellipsoid-level domain check, shared by the per-step gate and the endpoint gate.
function _domain_reason(E, X)
    if X isa LazySets.AbstractHyperrectangle
        _ellipsoid_in_box(E, X) && return nothing
        return _box_violation_reason(E, X)
    end

    if X isa UT.SetMinus
        included = UT.minus_included(X)
        included isa LazySets.AbstractHyperrectangle || return nothing
        _ellipsoid_in_box(E, included) || return _box_violation_reason(E, included)
        for hole in _each_member(UT.minus_hole(X))
            _provably_disjoint(E, hole) ||
                return "is not provably disjoint from a domain hole (obstacle)"
        end
        return nothing
    end

    return nothing
end

function state_domain_gate(rec::StepRecord, X)
    rec.status == :ok || return nothing
    reason = _domain_reason(rec.ellipsoid, X)
    reason === nothing && return nothing
    return "funnel ellipsoid at k=$(rec.k) " * reason
end

# ------------------------------------------------------------
# Endpoint gate. The chain's data endpoint — the backward terminal E_{K+1} or the
# forward entry E_1 — enters no StepRecord, so the per-step gates never see it;
# without this gate a certificate could end in (or start from) an ellipsoid that
# leaves the domain or crosses an obstacle.
# ------------------------------------------------------------

function endpoint_gate(E, X, r_min::Float64, check_state_domain::Bool, label::String)
    if r_min > 0.0
        Q = LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(E)))
        rmin = sqrt(max(0.0, LA.eigmin(Q)))
        rmin >= r_min || return "$label ellipsoid collapsed: min semi-axis " *
               "$(round(rmin; sigdigits = 3)) < r_min $(r_min)"
    end
    if check_state_domain
        reason = _domain_reason(E, X)
        reason === nothing || return "$label ellipsoid " * reason
    end
    return nothing
end

# ------------------------------------------------------------
# Gate D — terminal containment and the default terminal ellipsoid. The chain's
# guarantee ends in E_terminal; unless E_terminal ⊆ target_set the certificate does
# not certify the reach specification.
# ------------------------------------------------------------

# Largest axis-aligned ellipsoid centered at xN inside the target box (or inside the
# union member containing xN), shrunk by `shrink`. Returns (ellipsoid, nothing) or
# (nothing, reason).
function _default_terminal_ellipsoid(target_set, xN, shrink::Float64)
    # Prefer the box-centered inscription (maximal seed for the backward chain —
    # an endpoint-centered one is nearly minimal when the trajectory was truncated
    # at its first target hit, i.e. on the boundary). But the box-centered
    # ellipsoid excludes the box's corner shell: it is only usable when the
    # endpoint actually lies well inside it (margin 0.8), else the last backward
    # transition may be unreachable — fall back to the endpoint-centered
    # inscription. Ending trajectories centrally (deepest-hit truncation) is what
    # unlocks the big seed.
    for member in _each_member(target_set)
        member isa LazySets.AbstractHyperrectangle || continue
        v = Vector{Float64}(collect(xN))
        v ∈ member || continue
        c = Vector{Float64}(LazySets.center(member))
        r =
            shrink .*
            [LazySets.radius_hyperrectangle(member, i) for i in 1:LazySets.dim(member)]
        all(r .> 0) || continue
        if sum(((v .- c) ./ r) .^ 2) <= 0.8^2
            return LazySets.Ellipsoid(c, Matrix(LA.Diagonal(r .^ 2))), nothing
        end
        lo = LazySets.low(member)
        hi = LazySets.high(member)
        re = shrink .* min.(hi .- v, v .- lo)
        all(re .> 0) || continue
        return LazySets.Ellipsoid(v, Matrix(LA.Diagonal(re .^ 2))), nothing
    end
    return nothing,
    "no default terminal ellipsoid: the trajectory endpoint is not strictly inside " *
    "any hyperrectangle member of the target set (set options.terminal_shape " *
    "explicitly, or fix the trajectory)"
end

function terminal_containment(E_T, target_set)
    for member in _each_member(target_set)
        member isa LazySets.AbstractHyperrectangle || continue
        _ellipsoid_in_box(E_T, member) && return true
    end
    return false
end

# ------------------------------------------------------------
# Initial coverage (reported, never gating: the backward chain cannot control it).
# `initial_set ⊆ E_1` iff every vertex is inside (convexity); the margin is the max
# of (v−c)ᵀP(v−c) over the vertices — ≤ 1 means covered.
# ------------------------------------------------------------

function initial_coverage(initial_set, E1)
    E1 === nothing && return nothing
    verts = try
        LazySets.vertices_list(initial_set)
    catch
        return nothing
    end
    isempty(verts) && return nothing
    c = collect(Float64, LazySets.center(E1))
    P = Matrix{Float64}(UT.get_quadratic_form(E1))
    return maximum(LA.dot(v - c, P * (v - c)) for v in collect.(verts))
end

# ------------------------------------------------------------
# Gate application: convert a gate failure into an :infeasible record carrying the
# reason (the funnel and controller are dropped — they are not certified).
# ------------------------------------------------------------

function _gate_failure(rec::StepRecord, gate::Symbol, reason::String)
    return StepRecord(
        rec.k,
        :infeasible,
        nothing,
        nothing,
        nothing,
        merge(rec.summary, (; gate = gate, gate_reason = reason)),
    )
end

# Backward records carry the box-consistency evidence (required_X/U radii);
# forward boxes are known before solving, so their records carry none and get
# the tube-inflation guard instead.
function apply_gates(rec::StepRecord, ctx::ChainContext)
    rec.status == :ok || return rec

    reason = box_consistency_gate(rec)
    reason === nothing || return _gate_failure(rec, :box_consistency, reason)

    reason = collapse_gate(rec, ctx.options.r_min)
    reason === nothing || return _gate_failure(rec, :collapse, reason)

    if ctx.options.check_state_domain
        reason = state_domain_gate(rec, ctx.problem.system.X)
        reason === nothing || return _gate_failure(rec, :state_domain, reason)
    end

    return rec
end

function apply_gates(rec::StepRecord, problem, opts::ForwardOptions)
    rec.status == :ok || return rec

    reason = collapse_gate(rec, opts.r_min)
    reason === nothing || return _gate_failure(rec, :collapse, reason)

    if isfinite(opts.α_max) &&
       rec.summary.contraction isa Float64 &&
       !isnan(rec.summary.contraction) &&
       rec.summary.contraction > opts.α_max
        return _gate_failure(
            rec,
            :tube_inflation,
            "tube scale $(round(rec.summary.contraction; sigdigits = 3)) exceeds α_max",
        )
    end

    if opts.check_state_domain
        reason = state_domain_gate(rec, problem.system.X)
        reason === nothing || return _gate_failure(rec, :state_domain, reason)
    end

    return rec
end
