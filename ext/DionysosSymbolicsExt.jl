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

function _centered_box(center, radius)
    c = collect(Float64, center)
    r = collect(Float64, radius)
    @assert length(c) == length(r)
    return IA.interval.(c .- r, c .+ r)
end

function ST._getLipschitzConstants(J, xi, Xi_vals)
    L = zeros(length(xi))

    for (i, g) in enumerate(eachrow(J))
        Hg_s = Symbolics.jacobian(g, xi)
        nr, nc = size(Hg_s)

        mat = Matrix{Any}(undef, nr, nc)
        for r in 1:nr, c in 1:nc
            bf = Symbolics.build_function(
                Hg_s[r, c],
                xi...;
                expression = Val(false),
                nanmath = false,
            )

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

# Affine system + Lipschitz bounds around the linearization point, valid on the
# interval boxes X/U/W: symbolic Jacobians evaluated at the point; error
# constants bounded over the boxes with interval arithmetic.
function _affine_system_and_lipschitz(f, x, u, w, x̄, ū, w̄, X, U, W)
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

    function evalSym(z)
        bf = Symbolics.build_function(z, xi...; expression = Val(false), nanmath = false)

        fz = bf isa Tuple ? bf[1] : bf
        return Float64.(fz(x̄i...))
    end

    A = evalSym(Jx)
    B = evalSym(Ju)
    E = evalSym(Jw)
    f0 = vec(evalSym(f))

    c = f0 - A * x̄ - B * ū - E * w̄

    return MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W), L
end

# Automatic growth bound: trace the dynamics symbolically, bound each Jacobian
# entry over the state domain X with interval arithmetic (off-diagonal entries by
# magnitude, diagonal entries by their signed supremum, as required by the
# growth-bound ODE ṙ = L(u)·r), and return the bound as a function of the input.
function ST.compute_jacobian_bound(system::MS.ConstrainedBlackBoxControlContinuousSystem)
    nx = MS.statedim(system)
    X = MS.stateset(system)
    X === nothing && error(
        "compute_jacobian_bound requires the system to carry a bounded state set X; " *
        "provide `jacobian_bound` explicitly otherwise.",
    )

    xs = [Symbolics.variable(:x, i) for i in 1:nx]
    us = [Symbolics.variable(:u, i) for i in 1:MS.inputdim(system)]

    f = try
        collect(MS.mapping(system)(xs, us))
    catch err
        error(
            "compute_jacobian_bound: could not trace the dynamics symbolically " *
            "($(sprint(showerror, err))). Provide `jacobian_bound` explicitly.",
        )
    end

    J = Symbolics.jacobian(f, xs)

    LS = Dionysos.LazySets
    xint = [IA.interval(Float64(LS.low(X, i)), Float64(LS.high(X, i))) for i in 1:nx]

    entry_funs = Matrix{Any}(undef, nx, nx)
    for j in 1:nx, i in 1:nx
        bf = Symbolics.build_function(
            J[i, j],
            vcat(xs, us)...;
            expression = Val(false),
            nanmath = false,
        )
        entry_funs[i, j] = bf isa Tuple ? bf[1] : bf
    end

    _diag_bound(v::IA.Interval) = IA.sup(v)
    _diag_bound(v) = Float64(v)
    _offdiag_bound(v::IA.Interval) = IA.mag(v)
    _offdiag_bound(v) = abs(Float64(v))

    function jacobian_bound(u)
        args = vcat(xint, map(to_interval, collect(u)))
        M = zeros(nx, nx)
        for j in 1:nx, i in 1:nx
            v = entry_funs[i, j](args...)
            M[i, j] = i == j ? _diag_bound(v) : _offdiag_bound(v)
        end
        return ST.SMatrix{nx, nx}(M)
    end
    return jacobian_bound
end

function ST.build_affine_approximation(
    provider::ST.SymbolicAffineApproximationProvider,
    xbar,
    ubar,
    wbar = nothing;
    δx,
    δu,
)
    wbar === nothing && (wbar = zeros(length(provider.w)))

    Xbar = _centered_box(xbar, δx)
    Ubar = _centered_box(ubar, δu)
    Wbar = _centered_box(wbar, provider.ΔW)

    affineSys, L = _affine_system_and_lipschitz(
        provider.fsymbolic,
        provider.x,
        provider.u,
        provider.w,
        xbar,
        ubar,
        wbar,
        Xbar,
        Ubar,
        Wbar,
    )

    return ST.AffineApproximation(
        affineSys,
        L,
        provider.Uformat,
        provider.Wformat,
        (; δx = copy(δx), δu = copy(δu)),
    )
end

end # module
