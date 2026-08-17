# Step-level geometry helpers and diagnostics summaries. The radii helpers ARE
# load-bearing: `_required_linearization_box_radii` / `_box_contains_required_radii`
# back the box-consistency gate, `_ellipsoid_axis_radii` sizes the forward
# linearization box, and `_controller_image_axis_radii` backs the forward u-side
# consistency check. Only the `_step_summary` machinery is pure reporting.

function _ellipsoid_axis_radii(E)
    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    return sqrt.(max.(0.0, LA.diag(Q)))
end

# log√det(Q) = −log√det(P): larger is a bigger ellipsoid.
function _ellipsoid_logvolume(E)
    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    F = LA.cholesky(LA.Symmetric(Q))
    return sum(log, LA.diag(F.U))
end

function _controller_matrices(kappa::MS.AffineMap, nx::Int)
    return Matrix{Float64}(kappa.A), vec(Float64.(kappa.c))
end

function _controller_matrices(kappa::AbstractMatrix, nx::Int)
    K = Matrix{Float64}(kappa[:, 1:nx])
    b = vec(Float64.(kappa[:, nx + 1]))
    return K, b
end

function _controller_image_axis_radii(kappa, E)
    c = LazySets.center(E)
    K, b = _controller_matrices(kappa, length(c))

    Q = Matrix{Float64}(LazySets.shape_matrix(E))
    Σu = K * Q * K'
    uc = K * collect(Float64, c) + b
    ru = sqrt.(max.(0.0, LA.diag(LA.Symmetric(Σu))))

    return uc, ru
end

# Axis-aligned radii the linearization box must cover: the certified ellipsoid and
# its controller image, both measured from the linearization point.
function _required_linearization_box_radii(E_prev, kappa, xk, uk)
    rx = _ellipsoid_axis_radii(E_prev)
    uc, ru = _controller_image_axis_radii(kappa, E_prev)

    required_δx =
        abs.(collect(Float64, LazySets.center(E_prev)) .- collect(Float64, xk)) .+ rx
    required_δu = abs.(uc .- collect(Float64, uk)) .+ ru

    return required_δx, required_δu
end

function _box_contains_required_radii(δx, δu, required_δx, required_δu; atol::Float64)
    return all(required_δx .<= δx .+ atol) && all(required_δu .<= δu .+ atol)
end

function _candidate_diagnostics_empty()
    return (;
        candidate_scales = Float64[],
        candidate_logvolumes = Union{Nothing, Float64}[],
        candidate_statuses = Symbol[],
        candidate_Xbar_radii = Vector{Float64}[],
        candidate_Ubar_radii = Vector{Float64}[],
    )
end

function _step_summary(
    L,
    δx,
    δu,
    required_X_radius,
    required_U_radius,
    adaptive_box_iters,
    adaptive_box_status;
    selected_logvolume = nothing,
    selected_scale = nothing,
    selected_candidate_index = nothing,
    number_of_candidate_boxes = 0,
    candidate_scales = Float64[],
    candidate_logvolumes = Union{Nothing, Float64}[],
    candidate_statuses = Symbol[],
    candidate_Xbar_radii = Vector{Float64}[],
    candidate_Ubar_radii = Vector{Float64}[],
    provider_summary = nothing,
)
    return (;
        L,
        Xbar_radius = copy(δx),
        Ubar_radius = copy(δu),
        required_X_radius,
        required_U_radius,
        adaptive_box_iters,
        adaptive_box_status,
        selected_logvolume,
        selected_scale,
        selected_candidate_index,
        number_of_candidate_boxes,
        candidate_scales,
        candidate_logvolumes,
        candidate_statuses,
        candidate_Xbar_radii,
        candidate_Ubar_radii,
        provider_summary,
    )
end

function _append_candidate_diagnostic!(diag, scale, δx, δu, result)
    push!(diag.scales, Float64(scale))
    push!(diag.logvolumes, result.logvolume)
    push!(diag.statuses, result.status)
    push!(diag.Xbar_radii, copy(δx))
    push!(diag.Ubar_radii, copy(δu))
    return diag
end

function _candidate_diagnostics_tuple(diag)
    return (;
        candidate_scales = copy(diag.scales),
        candidate_logvolumes = copy(diag.logvolumes),
        candidate_statuses = copy(diag.statuses),
        candidate_Xbar_radii = copy(diag.Xbar_radii),
        candidate_Ubar_radii = copy(diag.Ubar_radii),
    )
end
