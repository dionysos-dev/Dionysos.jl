# ============================================================
# Main automatic level builder
# ============================================================

function build_levels_and_terminal_set(
    pclf::PCLF.PCLF,
    X,
    regions_to_avoid;
    tol::Float64 = 1e-3,
    max_levels::Int = 200,
    ΓX::Union{Nothing, Float64} = nothing,
    nb_levels::Union{Nothing, Int} = nothing,
)
    levels = build_levels_from_problem(
        pclf,
        X,
        regions_to_avoid;
        tol = tol,
        max_levels = max_levels,
        ΓX = ΓX,
        nb_levels = nb_levels,
    )

    D = compute_D_from_tau(pclf, levels[1])

    return levels, D
end

function build_levels_from_problem(
    pclf::PCLF.PCLF,
    X::Hyperrectangle,
    regions_to_avoid;
    tol::Float64 = 1e-3,
    max_levels::Int = 200,
    ΓX::Union{Nothing, Float64} = nothing,
    nb_levels::Union{Nothing, Int} = nothing,
)
    γ = pclf.JSRapprox
    0.0 < γ < 1.0 || error("Expected contraction rate γ in (0,1), got $γ.")
    max_levels >= 1 || error("Expected max_levels >= 1, got $max_levels.")
    tol >= 0.0 || error("Expected tol >= 0, got $tol.")

    if nb_levels !== nothing
        nb_levels >= 1 || error("Expected nb_levels >= 1, got $nb_levels.")
    end

    τX = if isnothing(ΓX)
        compute_tau_X(pclf, X; safety_factor = 1 + tol)
    else
        Float64(ΓX)
    end

    levels = Float64[]
    τ = Float64(τX)

    if nb_levels !== nothing
        for _ in 1:nb_levels
            push!(levels, τ)
            τ *= γ
        end
        reverse!(levels)
        return levels
    end

    τD = nothing
    for _ in 1:max_levels
        push!(levels, τ)

        if all_nodes_clear_regions(pclf, τ, regions_to_avoid; tol = tol)
            τD = τ
            break
        end

        τ *= γ
    end

    isnothing(τD) && error(
        "Could not find τD within max_levels=$max_levels. " *
        "Try increasing max_levels or relaxing the stopping criterion.",
    )

    reverse!(levels)
    return levels
end

# ============================================================
# τX computation
# ============================================================

function compute_tau_X(pclf::PCLF.PCLF, X::Hyperrectangle; safety_factor::Float64 = 1.05)
    safety_factor >= 1.0 || error("Expected safety_factor >= 1, got $safety_factor.")
    return safety_factor * gamma_cover_region_all_nodes(pclf, X)
end

function gamma_cover_region_all_nodes(pclf::PCLF.PCLF, X::Hyperrectangle)
    vals = Float64[]

    for piece in values(pclf.pieces)
        push!(vals, gamma_cover_set(piece, X))
    end

    isempty(vals) && error("PCLF has no pieces.")
    return maximum(vals)
end

"""
    gamma_cover_set(piece::PCLF.PolyhedralPiece, X::Hyperrectangle)

Return the smallest τ such that

    X ⊆ P_piece(τ)

for

    P_piece(τ) = {x : -τ w <= Gx <= τ w}.
"""
function gamma_cover_set(piece::PCLF.PolyhedralPiece, X::Hyperrectangle)
    G = piece.G
    w = piece.w

    vals = Float64[]
    for i in 1:size(G, 1)
        wi = Float64(w[i])
        wi > 0 || error("Expected strictly positive weight w[$i], got $wi.")

        gi = vec(Float64.(G[i, :]))
        worst = _support_abs_row_on_hyperrectangle(gi, X)
        push!(vals, worst / wi)
    end

    return isempty(vals) ? 0.0 : maximum(vals)
end

"""
    _support_abs_row_on_hyperrectangle(g, X)

Compute

    max_{x ∈ X} |g' x|

for a hyperrectangle `X`.
"""
function _support_abs_row_on_hyperrectangle(g::AbstractVector, X::Hyperrectangle)
    c = LazySets.center(X)
    r = radius_hyperrectangle(X)
    return abs(LA.dot(g, c)) + sum(abs.(g) .* r)
end

function radius_hyperrectangle(X::Hyperrectangle)
    return LazySets.radius_hyperrectangle(X)
end

# ============================================================
# Region intersection tests
# ============================================================

function all_nodes_clear_regions(
    pclf::PCLF.PCLF,
    τ::Real,
    regions_to_avoid;
    tol::Float64 = 0.0,
)
    for piece in values(pclf.pieces)
        for R in regions_to_avoid
            if piece_intersects_region_at_level(piece, τ, R; tol = tol)
                return false
            end
        end
    end
    return true
end

function piece_intersects_region_at_level(
    piece::PCLF.PolyhedralPiece,
    τ::Real,
    R;
    tol::Float64 = 0.0,
)
    tol >= 0.0 || error("Expected tol >= 0, got $tol.")

    τeff = max(Float64(τ) - tol, 0.0)
    Pτ = PCLF.get_sublevel_set(piece, τeff)
    I = set_intersection(Pτ, R)
    return is_nonempty_set(I)
end

# ============================================================
# Terminal-set construction
# ============================================================

function compute_D_from_tau(pclf::PCLF.PCLF, τD::Real)
    pieces = collect(values(pclf.pieces))
    isempty(pieces) && error("PCLF has no pieces.")

    sublevel_sets = [PCLF.get_sublevel_set(piece, τD) for piece in pieces]
    return SemiLinearSet(sublevel_sets)
end
