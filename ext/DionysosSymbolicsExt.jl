module DionysosSymbolicsExt

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

import Symbolics
import IntervalArithmetic as IA
import MathematicalSystems as MS
import LinearAlgebra as LA

function ST.interval_matrix_max_eig(mat::AbstractMatrix{<:IA.Interval})
    n, m = size(mat)
    @assert n == m "Matrix must be square"

    M = Array{Float64}(undef, n, n)
    for j in 1:n, i in 1:n
        a = mat[i, j]
        lo = IA.inf(a)
        hi = IA.sup(a)
        M[i, j] = max(abs(lo), abs(hi))
    end

    norm1 = maximum(sum(abs, M; dims = 1))[]
    normInf = maximum(sum(abs, M; dims = 2))[]

    return sqrt(norm1 * normInf)
end

to_interval(x) = x isa IA.Interval ? x : IA.interval(float(x), float(x))

function _bounds_to_intervals(lb, ub)
    @assert length(lb) == length(ub) "Lower/upper bounds must have the same dimension"
    return [IA.interval(float(lb[i]), float(ub[i])) for i in eachindex(lb)]
end

function _as_interval_vector(X)
    vals = collect(X)

    if all(v -> v isa IA.Interval, vals) || all(v -> v isa Real, vals)
        return map(to_interval, vals)
    end

    if length(vals) == 2
        lb, ub = vals[1], vals[2]

        if (lb isa Tuple || lb isa AbstractVector) &&
           (ub isa Tuple || ub isa AbstractVector)
            return _bounds_to_intervals(lb, ub)
        end
    end

    return error(
        "Cannot convert object of type $(typeof(X)) to a vector of intervals. " *
        "Pass either a vector of intervals, a vector of reals, or a box as [lb, ub].",
    )
end

function ST._getLipschitzConstants(J, xi, Xi_vals)
    L = zeros(length(xi))

    for (i, g) in enumerate(eachrow(J))
        Hg_s = Symbolics.jacobian(g, xi)
        nr, nc = size(Hg_s)

        mat = Matrix{Any}(undef, nr, nc)
        for r in 1:nr, c in 1:nc
            bf = Symbolics.build_function(Hg_s[r, c], xi...; expression = Val(false))
            f_rc = bf isa Tuple ? bf[1] : bf
            mat[r, c] = f_rc(Xi_vals...)
        end

        if any(x -> x isa IA.Interval, mat)
            matI = map(to_interval, mat)
            L[i] = ST.interval_matrix_max_eig(matI)
        else
            L[i] = maximum(abs, LA.eigen(Float64.(mat)).values)
        end
    end

    return L
end

function ST.buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)
    n = length(x)
    m = length(u)
    p = length(w)

    xi = vcat(x, u, w)
    x̄i = vcat(x̄, ū, w̄)

    Jx = Symbolics.jacobian(f, x)
    Ju = Symbolics.jacobian(f, u)
    Jw = Symbolics.jacobian(f, w)
    Jxi = hcat(Jx, Ju, Jw)

    Xi_vals = vcat(_as_interval_vector(X), _as_interval_vector(U), _as_interval_vector(W))
    L = ST._getLipschitzConstants(Jxi, xi, Xi_vals)

    sub_rules_x̄i = Dict(xi[i] => x̄i[i] for i in 1:(n + m + p))
    function evalSym(z)
        if z isa Number
            return Float64(Symbolics.value(Symbolics.substitute(z, sub_rules_x̄i)))
        else
            y = Symbolics.substitute.(z, Ref(sub_rules_x̄i))
            return Float64.(Symbolics.value.(y))
        end
    end

    A = evalSym(Jx)
    B = evalSym(Ju)
    E = evalSym(Jw)
    c = vec(evalSym(f) - A * x̄ - B * ū - E * w̄)

    return (MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W), L)
end

end
