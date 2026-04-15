struct EllipsoidalAffineApproximatedSystem
    dynamics::Dict{UT.Ellipsoid, MS.NoisyConstrainedAffineControlDiscreteSystem}
    L::Dict{UT.Ellipsoid, Float64}
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

function interval_matrix_max_eig(args...)
    return error(
        "interval_matrix_max_eig requires Symbolics.jl. " *
        "Load it with `using Symbolics` to enable this feature.",
    )
end

function _getLipschitzConstants(args...)
    return error(
        "_getLipschitzConstants requires Symbolics.jl. " *
        "Load it with `using Symbolics` to enable this feature.",
    )
end

function buildAffineApproximation(args...)
    return error(
        "buildAffineApproximation requires Symbolics.jl. " *
        "Load it with `using Symbolics` to enable this feature.",
    )
end
