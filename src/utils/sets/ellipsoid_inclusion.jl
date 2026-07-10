# Exact ellipsoid-in-ellipsoid inclusion test (an analytic kernel LazySets
# lacks), as the Dionysos-owned verb `is_included` — a `Base.issubset` method
# on two LazySets-owned types would be piracy. The math works in the
# quadratic-form matrices P = Q⁻¹, inverted once at entry.
function is_included(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
    P1 = get_quadratic_form(E1)
    P2 = get_quadratic_form(E2)
    c1 = LazySets.center(E1)
    c2 = LazySets.center(E2)

    e_min = eigmin(P1 - P2)
    if e_min < 0
        return false
    elseif c1 == c2
        return e_min >= 0
    elseif !(c1 ∈ E2)
        return false
    else
        L = cholesky(P2).L
        P = L \ P1 / L'
        c = L' * (c1 - c2)
        specDecomp = eigen(P)
        vals = specDecomp.values
        ct = specDecomp.vectors' * c

        g(β) = -β + sum((β * vals ./ (1 .- β * vals)) .* (ct .^ 2))
        dg(β) = -1 + sum((vals ./ (1 .- β * vals) .^ 2) .* (ct .^ 2))
        ddg(β) = 2 * sum((vals .^ 2 ./ (1 .- β * vals) .^ 3) .* (ct .^ 2))

        polPos(β) = -g(β)
        dpolPos(β) = -dg(β)
        ddpolPos(β) = -ddg(β)

        lb = 1 / min(vals...)
        ub = 1 - norm(ct)^2
        if ub < lb
            return false
        end
        (val, _) = derivative_bisection(
            polPos,
            dpolPos,
            ddpolPos;
            interval = [lb + 1e-15, ub],
            verbose = false,
            stopIfNegative = false,
        )
        gstar = -val
        return gstar >= -1
    end
end

function get_ℓ_ast_inclusion(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
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

    g(β) = -β + sum((β * vals ./ (1 .- β * vals)) .* (ct .^ 2))
    dg(β) = -1 + sum((vals ./ (1 .- β * vals) .^ 2) .* (ct .^ 2))
    ddg(β) = 2 * sum((vals .^ 2 ./ (1 .- β * vals) .^ 3) .* (ct .^ 2))

    polPos(β) = -g(β)
    dpolPos(β) = -dg(β)
    ddpolPos(β) = -ddg(β)

    lb = 1 / min(vals...)
    ub = norm(ct)^2 - 1
    if ub < lb
        ub = lb * 2
        while dpolPos(ub) < 0
            ub *= 2
        end
    end

    (val, βstar) = derivative_bisection(
        polPos,
        dpolPos,
        ddpolPos;
        interval = [lb + 1e-15, ub],
        verbose = false,
        stopIfNegative = false,
    )
    gstar = -val
    return -gstar, βstar
end

# return Enew, such that Enew is a scaled version of E1 and E2 ⊆_0 Enew
function scale_for_inclusion_contact_point(E1::LazySets.Ellipsoid, E2::LazySets.Ellipsoid)
    gstar, _ = get_ℓ_ast_inclusion(E2, E1)
    return get_sublevel_set(E1, gstar)
end
