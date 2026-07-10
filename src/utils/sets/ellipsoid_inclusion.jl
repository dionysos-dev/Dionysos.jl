
function Base.in(x::AbstractVector, elli::Ellipsoid)
    return (x - elli.c)'elli.P * (x - elli.c) ≤ 1
end

# Exact ellipsoid-in-ellipsoid inclusion test (an analytic kernel LazySets
# lacks), under the ecosystem verb: `E1 ⊆ E2`.
function Base.issubset(elli1::Ellipsoid, elli2::Ellipsoid)
    e_min = eigmin(elli1.P - elli2.P)
    if e_min < 0
        return false
    elseif elli1.c == elli2.c
        return e_min >= 0
    elseif !(elli1.c ∈ elli2)
        return false
    else
        L = cholesky(elli2.P).L
        P = L \ elli1.P / L'
        c = L' * (elli1.c - elli2.c)
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

function get_ℓ_ast_inclusion(elli1::Ellipsoid, elli2::Ellipsoid)
    L = cholesky(elli2.P).L
    P = L \ elli1.P / L'
    c = L' * (elli1.c - elli2.c)
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
function scale_for_inclusion_contact_point(E1::Ellipsoid, E2::Ellipsoid)
    gstar, _ = get_ℓ_ast_inclusion(E2, E1)
    return E1 * gstar
end
