function is_intersected(elli1::Ellipsoid, elli2::Ellipsoid)
    if elli1.c == elli2.c
        return true
    elseif (elli1.c ∈ elli2) || (elli2.c ∈ elli1)
        return true
    else
        L = cholesky(elli2.P).L
        P = L \ elli1.P / L'
        c = L' * (elli1.c - elli2.c)
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
        (val, _) =
            dbisection(polPos, dpolPos, ddpolPos; interval = [lb, ub], verbose = false) #, stopIfPositive=true)
        g_star = -val
        return g_star <= 1
    end
end

# return the smallest 2-norm of the elli1' where elli1' is the ellipsoid elli1 after 
# the change of variable transforming elli2 in B(0,1).

function get_ℓ_ast_intersect(elli1::Ellipsoid, elli2::Ellipsoid)
    L = cholesky(elli2.P).L
    P = L \ elli1.P / L'
    c = L' * (elli1.c - elli2.c)
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
    (val, βstar) = dbisection(
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
function scale_for_noninclusion_contact_point(E1::Ellipsoid, E2::Ellipsoid)
    if E1.c ∈ E2
        return nothing
    else
        gstar, _ = get_ℓ_ast_intersect(E2, E1)
        return E1 * gstar
    end
end
