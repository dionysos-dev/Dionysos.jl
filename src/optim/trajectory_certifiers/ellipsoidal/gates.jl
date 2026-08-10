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

function _ellipsoid_in_box(E, H::LazySets.AbstractHyperrectangle)
    c = collect(Float64, LazySets.center(E))
    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    cH = collect(Float64, LazySets.center(H))
    rH = [LazySets.radius_hyperrectangle(H, i) for i in 1:LazySets.dim(H)]
    return all(abs.(c .- cH) .+ sqrt.(max.(0.0, LA.diag(Q))) .<= rH .+ 1e-12)
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

function state_domain_gate(rec::StepRecord, X)
    rec.status == :ok || return nothing
    E = rec.ellipsoid

    if X isa LazySets.AbstractHyperrectangle
        _ellipsoid_in_box(E, X) && return nothing
        return "funnel ellipsoid at k=$(rec.k) is not contained in the state domain"
    end

    if X isa UT.SetMinus
        included = UT.minus_included(X)
        included isa LazySets.AbstractHyperrectangle || return nothing
        _ellipsoid_in_box(E, included) ||
            return "funnel ellipsoid at k=$(rec.k) is not contained in the state domain"
        for hole in _each_member(UT.minus_hole(X))
            _provably_disjoint(E, hole) ||
                return "funnel ellipsoid at k=$(rec.k) is not provably disjoint " *
                       "from a domain hole (obstacle)"
        end
        return nothing
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
    for member in _each_member(target_set)
        member isa LazySets.AbstractHyperrectangle || continue
        collect(xN) ∈ member || continue
        lo = LazySets.low(member)
        hi = LazySets.high(member)
        r = shrink .* min.(hi .- collect(xN), collect(xN) .- lo)
        all(r .> 0) || continue
        return LazySets.Ellipsoid(collect(float.(xN)), Matrix(LA.Diagonal(r .^ 2))), nothing
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
