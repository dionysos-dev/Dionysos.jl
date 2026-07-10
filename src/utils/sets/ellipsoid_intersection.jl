# Exact ellipsoid-ellipsoid disjointness (an analytic kernel LazySets lacks),
# as the Dionysos-owned verb `is_disjoint` — a `Base.isdisjoint` method on two
# LazySets-owned types would be piracy. The math works in the quadratic-form
# matrices P = Q⁻¹, inverted once at entry.
is_disjoint(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid) = !_intersects(E1, E2)

function _intersects(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
    c1 = LazySets.center(E1)
    c2 = LazySets.center(E2)
    if c1 == c2
        return true
    elseif (c1 ∈ E2) || (c2 ∈ E1)
        return true
    else
        P1 = get_quadratic_form(E1)
        P2 = get_quadratic_form(E2)
        L = cholesky(P2).L
        P = L \ P1 / L'
        c = L' * (c1 - c2)
        specDecomp = eigen(P)
        vals = specDecomp.values
        ct = specDecomp.vectors' * c

        g(β) = -β + sum((β * vals ./ (1 .+ β * vals)) .* (ct .^ 2))
        dg(β) = -1 + sum((vals ./ (1 .+ β * vals) .^ 2) .* (ct .^ 2))
        ddg(β) = -2 * sum((vals .^ 2 ./ (1 .+ β * vals) .^ 3) .* (ct .^ 2))
        polPos(β) = -g(β)
        dpolPos(β) = -dg(β)
        ddpolPos(β) = -ddg(β)
        lb = 0.0
        ub = norm(ct)^2 - 1
        if (ub < lb)
            return true
        end
        (val, _) = derivative_bisection(
            polPos,
            dpolPos,
            ddpolPos;
            interval = [lb, ub],
            verbose = false,
        )
        g_star = -val
        return g_star <= 1
    end
end

# return the smallest 2-norm of the elli1' where elli1' is the ellipsoid elli1 after
# the change of variable transforming elli2 in B(0,1).
function get_ℓ_ast_intersect(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
    P1 = get_quadratic_form(E1)
    P2 = get_quadratic_form(E2)
    c1 = LazySets.center(E1)
    c2 = LazySets.center(E2)

    L = cholesky(P2).L
    P = L \ P1 / L'
    c = L' * (c1 - c2)
    specDecomp = eigen(P)
    vals = specDecomp.values
    ct = specDecomp.vectors' * c

    g(β) = -β + sum((β * vals ./ (1 .+ β * vals)) .* (ct .^ 2))
    dg(β) = -1 + sum((vals ./ (1 .+ β * vals) .^ 2) .* (ct .^ 2))
    ddg(β) = -2 * sum((vals .^ 2 ./ (1 .+ β * vals) .^ 3) .* (ct .^ 2))
    polPos(β) = -g(β)
    dpolPos(β) = -dg(β)
    ddpolPos(β) = -ddg(β)

    lb = 0.0
    ub = norm(ct)^2 - 1
    if ub < lb
        ub = 10e-2
        while dpolPos(ub) < 0
            ub *= 2
        end
    end
    (val, βstar) = derivative_bisection(
        polPos,
        dpolPos,
        ddpolPos;
        interval = [lb, ub],
        verbose = false,
        stopIfNegative = false,
    )
    gstar = -val
    return gstar, βstar
end

# return Enew, such that Enew is a scaled version of E1 and E2 ∩_0 Enew
# return nothing if it is impossible, i.e. if E1.c ∈ E2
function scale_for_noninclusion_contact_point(
    E1::LazySets.Ellipsoid,
    E2::LazySets.Ellipsoid,
)
    if LazySets.center(E1) ∈ E2
        return nothing
    else
        gstar, _ = get_ℓ_ast_intersect(E2, E1)
        return get_sublevel_set(E1, gstar)
    end
end
