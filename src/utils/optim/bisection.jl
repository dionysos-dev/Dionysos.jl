
function bisection(
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
