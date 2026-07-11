# Helpers for `LazySets.Ellipsoid` — E(c, Q) = {x : (x−c)ᵀ Q⁻¹ (x−c) ≤ 1} with
# shape matrix `Q`. Dionysos's synthesis math (S-procedure LMIs, exact
# inclusion/intersection kernels) is written in the quadratic-form matrix
# P = Q⁻¹; `get_quadratic_form` bridges, and the analytic kernels invert once
# at entry.

_symmetrize(M) = (M + M') / 2

"""
    get_quadratic_form(E::LazySets.Ellipsoid)

Quadratic-form matrix `P = Q⁻¹` of the ellipsoid, i.e. the matrix such that
`E = {x : (x − c)ᵀ P (x − c) ≤ 1}`.
"""
get_quadratic_form(E::LazySets.Ellipsoid) = _symmetrize(inv(LazySets.shape_matrix(E)))

center_distance(S1, S2) = norm(LazySets.center(S1) - LazySets.center(S2))
center_distance(S, x::AbstractVector) = norm(LazySets.center(S) - x)

# Volume of {x : (x-c)'Q⁻¹(x-c) ≤ 1} = π^(N/2) / Γ(N/2 + 1) · √det(Q)
# (the volume of the unit N-ball scaled by the ellipsoid's semi-axes).
function get_volume(E::LazySets.Ellipsoid)
    N = LazySets.dim(E)
    return pi^(N / 2) / SpecialFunctions.gamma(N / 2 + 1) *
           sqrt(det(LazySets.shape_matrix(E)))
end

"Sublevel-set scaling `{x : (x−c)ᵀP(x−c) ≤ α}`, i.e. `Q ← α·Q`."
function get_sublevel_set(E::LazySets.Ellipsoid, α)
    return LazySets.Ellipsoid(
        LazySets.center(E),
        _symmetrize(α * LazySets.shape_matrix(E));
        check_posdef = false,
    )
end

"Semi-axis lengths `√λᵢ(Q)` of the ellipsoid, longest first."
function get_length_semiaxis(E::LazySets.Ellipsoid)
    return sqrt.(
        sort(eigen(Symmetric(Matrix(LazySets.shape_matrix(E)))).values; rev = true),
    )
end

include("ellipsoid_inclusion.jl")
include("ellipsoid_intersection.jl")

# compress E1 if E1∩E2≠∅
# return nothing if impossible
function compress_if_intersection(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
    if !is_disjoint(E1, E2)
        return scale_for_noninclusion_contact_point(E1, E2)
    else
        return E1
    end
end

# For a non-ellipsoidal obstacle there is no analytic shrink: keep E1 when it
# already avoids the obstacle, reject it otherwise.
compress_if_intersection(E1::LazySets.Ellipsoid, O) = is_disjoint(E1, O) ? E1 : nothing
