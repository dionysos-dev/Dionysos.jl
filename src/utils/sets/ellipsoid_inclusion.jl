
function Base.in(x::AbstractVecOrMat, elli::Ellipsoid)
    return (x - elli.c)'elli.P * (x - elli.c) ≤ 1
end

function Base.:∉(elli1::Ellipsoid, elli2::Ellipsoid)
    return !(elli1 ∈ elli2)
end

function Base.in(elli1::Ellipsoid, elli2::Ellipsoid; eps = 1e-10)
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
        (val, _) = dbisection(
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

function get_ℓ_ast_inclusion(elli1::Ellipsoid, elli2::Ellipsoid; eps = 1e-10)
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

    (val, βstar) = dbisection(
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

# minimise a convex function
function dbisection(
    f::Function,
    df::Function,
    ddf::Function;
    interval = [0, 1],
    δ = 1e-8,
    verbose = false,
    stopIfNegative = false,
)
    x = zeros(2)
    fval = zeros(2)
    dfval = zeros(2)
    L = 0

    x[1] = interval[1]
    x[2] = interval[2]

    fval[1] = f(x[1])
    fval[2] = f(x[2])
    dfval[1] = df(x[1])
    dfval[2] = df(x[2])
    L = ddf(x[1])
    if dfval[2] < 0
        return (fval[2], x[2])
    end

    fβ, β = let β = (x[1] + x[2]) / 2, fβ = f(β)
        dfβ = df(β)
        if (verbose)
            println(string(dfval[1]) * " " * string(dfval[2]) * "    " * string(fβ))
        end
        while abs(x[1] - x[2]) > δ &&
            (!stopIfNegative || (fβ > 0 && 2 * fβ < -L * (x[2] - x[1])^2))
            if (verbose)
                println(string(x[1]) * "\t" * string(x[2]) * "\t" * string(fβ))
                println(string(dfval[1]) * "\t" * string(dfval[2]) * "\t" * string(fβ))
            end
            if (dfβ < 0)
                x[1] = β
                dfval[1] = dfβ
                L = ddf(β)
            else
                x[2] = β
                dfval[2] = dfβ
            end

            β = (x[1] + x[2]) / 2
            fβ = f(β)
            dfβ = df(β)
        end
        (fβ, β)
    end
    return (fβ, β)
end
