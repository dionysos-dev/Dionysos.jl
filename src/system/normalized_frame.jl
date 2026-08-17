# Geometry of the globally normalized certification frame z = x ./ t
# (plan.md §4.3). The provider half lives in the Symbolics extension
# (`normalized_symbolic_provider`, which builds the scaled dynamics once so the
# Hessian bounds are exact in the working frame); this file is the set/trajectory/
# controller half every consumer used to hand-roll.

import LazySets

"""
    normalize_box(H, t) -> LazySets.Hyperrectangle

The axis-aligned box `H` in the normalized frame `z = x ./ t`.
"""
function normalize_box(H::LazySets.AbstractHyperrectangle, t::AbstractVector)
    return LazySets.Hyperrectangle(;
        low = collect(LazySets.low(H) ./ t),
        high = collect(LazySets.high(H) ./ t),
    )
end

"""
    normalize_trajectory(traj::Trajectory, t) -> Trajectory

The trajectory with every state scaled to the normalized frame `z = x ./ t`
(inputs are frame-invariant).
"""
function normalize_trajectory(traj::Trajectory, t::AbstractVector)
    n = length(t)
    return Trajectory(
        [SVector{n}(collect(x) ./ t) for x in states(traj)];
        inputs = collect(inputs(traj)),
    )
end

"""
    denormalize_ellipsoid(E, t) -> LazySets.Ellipsoid

An ellipsoid certified in the normalized frame, mapped back to the physical
frame: `c_x = t .* c_z`, `Q_x = D·Q_z·D` with `D = Diagonal(t)`.
"""
function denormalize_ellipsoid(E::LazySets.Ellipsoid, t::AbstractVector)
    D = LA.Diagonal(collect(Float64, t))
    return LazySets.Ellipsoid(
        t .* collect(LazySets.center(E)),
        Matrix(D * Matrix(LazySets.shape_matrix(E)) * D);
        check_posdef = false,
    )
end

"""
    denormalize_controller(κ::MathematicalSystems.AffineMap, t) -> AffineMap

A feedback `u = K_z·z + b` synthesized in the normalized frame, as the physical
map `u = K_z·D⁻¹·x + b` (inputs are frame-invariant).
"""
function denormalize_controller(κ::MS.AffineMap, t::AbstractVector)
    D = LA.Diagonal(collect(Float64, t))
    return MS.AffineMap(Matrix(Matrix(κ.A) / D), collect(κ.c))
end

"""
    denormalize_funnel(ctrl::FunnelController, t) -> FunnelController

A whole certified funnel controller mapped back to the physical frame — the
one-call version of denormalizing every `κ_k` and every `E_k`.
"""
function denormalize_funnel(ctrl::FunnelController, t::AbstractVector)
    return FunnelController(
        [denormalize_controller(κ, t) for κ in ctrl.kappas],
        [denormalize_ellipsoid(E, t) for E in ctrl.ellipsoids],
    )
end
