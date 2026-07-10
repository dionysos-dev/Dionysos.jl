# ============================================================
# Scalar optimization on an interval.
#
# Three small routines for minimizing a scalar function `f : [a, b] → ℝ`,
# differing in the derivative information they exploit:
#
#   golden_section_search  derivative-free (unimodal `f`)
#   newton_method          uses `f, f′, f″`
#   derivative_bisection   bisects on `f′` to locate the min of a convex `f`
#
# All return `(f(x⋆), x⋆)`. `stopIfNegative` enables the domain-specific early
# exit used by the ellipsoid inclusion/intersection tests (bail out as soon as
# the sign of `f` settles the query, before full convergence).
# ============================================================

"""
    golden_section_search(f; interval, δ, verbose, stopIfNegative) -> (f(x⋆), x⋆)

Golden-section search for the minimum of a unimodal `f` on `interval`. Derivative
free; brackets the minimum by keeping the golden ratio between four probe points.
"""
function golden_section_search(
    f::Function;
    interval = [0, 1],
    δ = 1e-8,
    verbose = false,
    stopIfNegative = false,
)
    phi = (1 + sqrt(5)) / 2
    x = zeros(4)
    fval = zeros(4)
    x[1] = interval[1]
    x[4] = interval[2]

    fval[1] = f(x[1])
    fval[4] = f(x[4])

    x[2] = x[4] - (x[4] - x[1]) / phi
    x[3] = x[1] + (x[4] - x[1]) / phi
    fval[2] = f(x[2])
    fval[3] = f(x[3])
    k = 0
    while abs(x[1] - x[4]) > δ && (!stopIfNegative || all(i -> i >= 0, fval))
        k = k + 1

        if (fval[2] < fval[3])
            x[4] = x[3]
            fval[4] = fval[3]

            x[3] = x[2]
            fval[3] = fval[2]

            x[2] = x[4] - (x[4] - x[1]) / phi
            fval[2] = f(x[2])
        else
            x[1] = x[2]
            fval[1] = fval[2]

            x[2] = x[3]
            fval[2] = fval[3]

            x[3] = x[1] + (x[4] - x[1]) / phi
            fval[3] = f(x[3])
        end
        if (verbose)
            println(min(fval...))
        end
    end
    (_, im) = findmin(fval)
    return (fval[im], x[im])
end

_projection(x, I) = min(max(x, I[1]), I[2])
_is_out(x, I) = x <= I[1] || x >= I[2]

"""
    newton_method(f, df, ddf; interval, x0, ϵ, verbose, stopIfNegative) -> (f(x⋆), x⋆)

Newton's method for a stationary point of `f`, iterating `x -= f′(x)/f″(x)` until
`|f′(x)| ≤ ϵ` or `x` leaves `interval`; the result is projected back onto `interval`.
"""
function newton_method(
    f::Function,
    df::Function,
    ddf::Function;
    interval = [0, 1],
    x0 = (interval[1] + interval[2]) / 2.0,
    ϵ = 1e-8,
    verbose = false,
    stopIfNegative = false,
)
    x = x0
    while (abs(df(x)) > ϵ && (!stopIfNegative || f(x) ≥ 0) && !_is_out(x, interval))
        x = x - (df(x) / ddf(x))
        if (verbose)
            println(f(x))
        end
    end
    x = _projection(x, interval)
    return (f(x), x)
end

"""
    derivative_bisection(f, df, ddf; interval, δ, verbose, stopIfNegative) -> (f(x⋆), x⋆)

Minimize a convex `f` by bisecting on the sign of `f′`: the sub-interval is halved
toward where the derivative changes sign. `ddf` supplies the local curvature bound
used by the `stopIfNegative` early exit.
"""
function derivative_bisection(
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
