function projection(x, I)
    return min(max(x, I[1]), I[2])
end

function is_out(x, I)
    return x <= I[1] || x >= I[2]
end

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
    while (abs(df(x)) > ϵ && (!stopIfNegative || f(x) ≥ 0) && !is_out(x, interval))
        x = x - (df(x) / ddf(x))
        if (verbose)
            println(f(x))
        end
    end
    x = projection(x, interval)
    return (f(x), x)
end
