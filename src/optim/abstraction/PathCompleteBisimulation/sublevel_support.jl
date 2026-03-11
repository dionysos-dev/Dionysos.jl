# ============================================================
# Automatic level selection
# ============================================================

"""
    _support_abs_row_on_hyperrectangle(g, X)

Compute max_{x in X} |g'x| for a hyperrectangle X.
"""
function _support_abs_row_on_hyperrectangle(g::AbstractVector, X::Hyperrectangle)
    c = center(X)
    r = radius_hyperrectangle(X)
    return abs(LA.dot(g, c)) + sum(abs.(g) .* r)
end

"""
    gamma_cover_set(piece::PCLF.PolyhedralPiece, X::Hyperrectangle)

Return the smallest Γ such that X ⊆ P_piece(Γ).
For P_piece(Γ) = {x : -Γ w <= Gx <= Γ w}.
"""
function gamma_cover_set(piece::PCLF.PolyhedralPiece, X::Hyperrectangle)
    G = piece.G
    w = piece.w

    vals = Float64[]
    for i in 1:size(G, 1)
        gi = vec(G[i, :])
        worst = _support_abs_row_on_hyperrectangle(gi, X)
        push!(vals, worst / w[i])
    end
    return maximum(vals)
end

"""
    gamma_cover_region_all_nodes(pclf, X)

Return Γ_X such that X ⊆ P_s(Γ_X) for all nodes s.
Currently implemented for polyhedral pieces and Hyperrectangle regions.
"""
function gamma_cover_region_all_nodes(pclf::PCLF.PCLF, X::Hyperrectangle)
    vals = Float64[]
    for piece in values(pclf.pieces)
        if piece isa PCLF.PolyhedralPiece
            push!(vals, gamma_cover_set(piece, X))
        else
            error("Automatic Γ computation currently implemented only for PolyhedralPiece.")
        end
    end
    return maximum(vals)
end

"""
    gamma_cover_terminal_all_nodes(pclf, D)

Return Γ_D such that D ⊆ P_s(Γ_D) for all nodes s.
"""
function gamma_cover_terminal_all_nodes(pclf::PCLF.PCLF, D::Hyperrectangle)
    vals = Float64[]
    for piece in values(pclf.pieces)
        if piece isa PCLF.PolyhedralPiece
            push!(vals, gamma_cover_set(piece, D))
        else
            error("Automatic Γ computation currently implemented only for PolyhedralPiece.")
        end
    end
    return maximum(vals)
end

"""
    build_levels_from_problem(pclf, X, D; num_levels = nothing, safety_factor = 1.05)

Automatically build the Lyapunov levels Γ if the user did not provide them.

Returns `(Γ, τ, ΓD, ΓX)` where:
- `ΓX` covers the full region X,
- `ΓD` covers the terminal region D,
- `τ` is the ratio between successive levels.

If `num_levels === nothing`, choose the smallest N such that
ΓD * (1/γ)^(N-1) >= ΓX, where γ = pclf.JSRapprox.
"""
function build_levels_from_problem(
    pclf::PCLF.PCLF,
    X::Hyperrectangle,
    D::Hyperrectangle;
    num_levels::Union{Nothing,Int} = nothing,
    safety_factor::Float64 = 1.05,
)
    γ = pclf.JSRapprox
    0.0 < γ < 1.0 || error("Expected contraction rate γ in (0,1), got $γ.")

    ΓX = safety_factor * gamma_cover_region_all_nodes(pclf, X)
    ΓD = safety_factor * gamma_cover_terminal_all_nodes(pclf, D)

    ΓD > 0 || error("Computed ΓD <= 0, cannot build automatic levels.")
    ΓX >= ΓD || error("Expected ΓX >= ΓD, got ΓX=$ΓX and ΓD=$ΓD.")

    # Default ratio induced by contraction
    τ_default = 1 / γ

    if isnothing(num_levels)
        if isapprox(ΓX, ΓD; atol=1e-14, rtol=1e-12)
            num_levels = 1
        else
            num_levels = Int(ceil(log(ΓX / ΓD) / log(τ_default))) + 1
        end
    end

    num_levels >= 1 || error("num_levels must be >= 1.")

    if num_levels == 1
        Γ = [ΓX]
        τ = 1.0
    else
        # exact interpolation between ΓD and ΓX
        τ = (ΓX / ΓD)^(1 / (num_levels - 1))
        Γ = [ΓD * τ^(i - 1) for i in 1:num_levels]
    end

    return Float64.(Γ), τ, ΓD, ΓX
end