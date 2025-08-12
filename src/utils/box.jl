function sample(X::IntervalBox)
    return Vector(map(x -> x.lo + (x.hi - x.lo) * rand(), X.v))
end

@recipe function f(X::IntervalBox; dims = [1, 2])
    xmin = Vector(map(x -> x.lo, X.v))
    xmax = Vector(map(x -> x.hi, X.v))
    dims := dims
    return HyperRectangle(xmin, xmax)
end
