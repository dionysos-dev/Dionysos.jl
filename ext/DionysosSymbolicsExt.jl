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

    # Worst-case absolute matrix |A|
    M = Array{Float64}(undef, n, n)
    for j in 1:n, i in 1:n
        a = mat[i, j]
        lo = IA.inf(a)
        hi = IA.sup(a)
        M[i, j] = max(abs(lo), abs(hi))
    end

    # 1-norm (max column sum)
    norm1 = maximum(sum(abs, M; dims = 1))[]
    # ∞-norm (max row sum)
    normInf = maximum(sum(abs, M; dims = 2))[]

    return sqrt(norm1 * normInf)
end

to_interval(x) = x isa IA.Interval ? x : IA.interval(float(x), float(x))
_to_interval_box(X) = X isa IA.IntervalBox ? X : IA.IntervalBox(X)

function ST._getLipschitzConstants(J, xi, Xi_vals)
    L = zeros(length(xi))

    for (i, g) in enumerate(eachrow(J))
        Hg_s = Symbolics.jacobian(g, xi)
        nr, nc = size(Hg_s)

        mat = Matrix{Any}(undef, nr, nc)
        for r in 1:nr, c in 1:nc
            f_rc =
                Symbolics.build_function(Hg_s[r, c], collect(xi); expression = Val(false))
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
    X = _to_interval_box(X)
    U = _to_interval_box(U)
    W = _to_interval_box(W)

    n = Base.length(x)
    m = Base.length(u)
    p = Base.length(w)

    xi = vcat(x, u, w)
    x̄i = vcat(x̄, ū, w̄)

    Jx = Symbolics.jacobian(f, x)
    Ju = Symbolics.jacobian(f, u)
    Jw = Symbolics.jacobian(f, w)
    Jxi = hcat(Jx, Ju, Jw)

    Xi_vals = vcat(collect(X), collect(U), collect(W))
    L = ST._getLipschitzConstants(Jxi, xi, Xi_vals)

    sub_rules_x̄i = Dict(xi[i] => x̄i[i] for i in 1:(n + m + p))
    function evalSym(x)
        if x isa Number
            return Float64(Symbolics.value(Symbolics.substitute(x, sub_rules_x̄i)))
        else
            y = Symbolics.substitute.(x, Ref(sub_rules_x̄i))
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
