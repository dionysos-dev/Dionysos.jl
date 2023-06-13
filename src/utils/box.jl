using IntervalArithmetic

function plot_box!(X::IntervalBox; dims = [1, 2], color = :red, opacity = 1.0)
    xmin = Vector(map(x -> x.lo, X.v))
    xmax = Vector(map(x -> x.hi, X.v))
    c = (xmin + xmax) ./ 2
    h = xmax - xmin
    return plot!(
        rectangle(c[dims], h ./ 2);
        opacity = opacity,
        color = color,
        legend = false,
    )
end

function sample_box(X::IntervalBox)
    return Vector(map(x -> x.lo + (x.hi - x.lo) * rand(), X.v))
end
