import IntervalLinearAlgebra as ILA

struct EllipsoidalAffineApproximatedSystem{}
    dynamics::Dict{UT.Ellipsoid, MS.NoisyConstrainedAffineControlDiscreteSystem}
    L::Dict{UT.Ellipsoid, Float64} # smoothness constant to bound error
end

# --------------------------
# Helperss
# --------------------------
to_interval(x) = x isa IA.Interval ? x : IA.interval(float(x), float(x))

# Force a matrix into a concretely typed interval matrix (avoids Matrix{Real}).
function as_interval_matrix(A)
    AI = map(to_interval, A)
    return Matrix{IA.Interval{Float64}}(AI)
end

# Conservative symmetric-ification (useful for Hessians).
symmetrize(A) = (A .+ transpose(A)) ./ 2

# Max absolute value of an interval.
absmax(a::IA.Interval) = max(abs(IA.inf(a)), abs(IA.sup(a)))

# Complex interval magnitude upper bound
# For z = a + i b with a,b intervals, |z| <= sqrt( max|a|^2 + max|b|^2 )
function absmax(z::Complex{<:IA.Interval})
    ar = absmax(real(z))
    ai = absmax(imag(z))
    return sqrt(ar^2 + ai^2)
end

# If eigenbox ever returns a plain complex number
absmax(z::Complex{<:Real}) = abs(z)

# --------------------------
# Method 1: Always-safe bound (very conservative)
# --------------------------

function interval_matrix_max_eig(mat::AbstractMatrix{<:IA.Interval})
    n, m = size(mat)
    @assert n == m "Matrix must be square"

    M = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        a = mat[i, j]
        M[i, j] = max(abs(IA.inf(a)), abs(IA.sup(a)))
    end

    norm1 = maximum(sum(abs, M; dims = 1))[]
    normInf = maximum(sum(abs, M; dims = 2))[]

    return sqrt(norm1 * normInf)
end

# --------------------------
# Method 2: Tighter bound using IntervalLinearAlgebra (when it works)
# --------------------------

function interval_matrix_eigenbox_bound(mat::AbstractMatrix)
    HI = Matrix{IA.Interval{Float64}}(
        map(x -> x isa IA.Interval ? x : IA.interval(float(x), float(x)), mat),
    )
    HI = (HI .+ transpose(HI)) ./ 2  # encourage real eigenvalue enclosure

    λ = ILA.eigenbox(HI)             # may be intervals OR complex intervals
    return maximum(absmax, λ)        # absmax now supports both
end

# --------------------------
# Unified Lipschitz constant computation
# --------------------------

function _getLipschitzConstants(J, xi, rules; use_eigenbox::Bool = true)
    L = zeros(Float64, length(xi))

    for (i, g) in enumerate(eachrow(J))
        Hg_s = Symbolics.jacobian(g, xi)          # Hessian (symbolic)
        Hg = Symbolics.substitute(Hg_s, rules)  # substitute intervals
        mat = Symbolics.value.(Hg)

        if any(x -> x isa IA.Interval, mat)
            # interval / mixed case
            HI = as_interval_matrix(mat)

            if use_eigenbox
                L[i] = interval_matrix_eigenbox_bound(HI)
            else
                L[i] = interval_matrix_max_eig(HI)
            end
        else
            # pure real case
            vals = LA.eigen(Matrix{Float64}(mat)).values
            L[i] = maximum(abs.(vals))
        end
    end

    return L
end

function buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)
    n = Base.length(x)
    m = Base.length(u)
    p = Base.length(w)
    xi = vcat(x, u, w)
    x̄i = vcat(x̄, ū, w̄)
    Xi = vcat(collect(X), collect(U), collect(W))
    sub_rules_Xi = Dict(xi[i] => Xi[i] for i in 1:(n + m + p))

    Jx = Symbolics.jacobian(f, x)
    Ju = Symbolics.jacobian(f, u)
    Jw = Symbolics.jacobian(f, w)
    Jxi = hcat(Jx, Ju, Jw)

    L = _getLipschitzConstants(Jxi, xi, sub_rules_Xi)

    sub_rules_x̄i = Dict(xi[i] => x̄i[i] for i in 1:(n + m + p))
    function evalSym(x)
        # When x is a Vector{SymbolicUtils.BasicSymbolic{Real}}, 
        # one needs to substitute each element of the vector
        if eltype(x) <: Symbolics.SymbolicUtils.BasicSymbolic{Real}
            return [Symbolics.substitute(elem, sub_rules_x̄i) for elem in x]
        end
        return Float64.(Symbolics.value.(Symbolics.substitute(x, sub_rules_x̄i)))
    end

    A = evalSym(Jx)
    B = evalSym(Ju)
    E = evalSym(Jw)
    c = vec(evalSym(f) - A * x̄ - B * ū - E * w̄)
    return (MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W), L)
end

struct AffineApproximationDiscreteSystem
    constrainedAffineSys::MS.NoisyConstrainedAffineControlDiscreteSystem
    L::Any
    f_eval::Any
    function AffineApproximationDiscreteSystem(sys, L)
        f_eval_fun(x, u, w) = sys.A * x + sys.B * u + sys.D * w + sys.c
        return new(sys, L, f_eval_fun)
    end
end
function AffineApproximationDiscreteSystem(A, B, c, E, X, U, W, L)
    contSys = MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W)
    return AffineApproximationDiscreteSystem(contSys, L)
end

struct SymbolicSystem
    fsymbolicT::Any
    fsymbolic::Any
    Ts::Any
    nx::Any
    nu::Any
    nw::Any
    x::Any
    u::Any
    w::Any
    ΔX::Any
    ΔU::Any
    ΔW::Any
    X::Any
    U::Any
    W::Any
    obstacles::Any
    f_eval::Any
    f_backward_eval::Any
    Uformat::Any
    Wformat::Any
end
