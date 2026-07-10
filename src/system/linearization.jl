struct EllipsoidalAffineApproximatedSystem
    dynamics::Dict{LazySets.Ellipsoid, MS.NoisyConstrainedAffineControlDiscreteSystem}
    L::Dict{LazySets.Ellipsoid, Float64}
end

struct AffineApproximationDiscreteSystem{
    S <: MS.NoisyConstrainedAffineControlDiscreteSystem,
    LT,
    F,
}
    constrainedAffineSys::S
    L::LT
    f_eval::F
    function AffineApproximationDiscreteSystem(
        sys::MS.NoisyConstrainedAffineControlDiscreteSystem,
        L,
    )
        f_eval_fun(x, u, w) = sys.A * x + sys.B * u + sys.D * w + sys.c
        return new{typeof(sys), typeof(L), typeof(f_eval_fun)}(sys, L, f_eval_fun)
    end
end

function AffineApproximationDiscreteSystem(A, B, c, E, X, U, W, L)
    contSys = MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W)
    return AffineApproximationDiscreteSystem(contSys, L)
end

# Bundle built by the Symbolics extension (`buildAffineApproximation`): symbolic
# dynamics + local approximation domains + evaluators. One type parameter per
# field keeps it type-stable despite the heterogeneous contents.
struct SymbolicSystem{
    FT,
    FS,
    TS,
    NX,
    NU,
    NW,
    XT,
    UTT,
    WT,
    DX,
    DU,
    DW,
    XS,
    US,
    WS,
    OB,
    FE,
    FB,
    UF,
    WF,
}
    fsymbolicT::FT
    fsymbolic::FS
    Ts::TS
    nx::NX
    nu::NU
    nw::NW
    x::XT
    u::UTT
    w::WT
    ΔX::DX
    ΔU::DU
    ΔW::DW
    X::XS
    U::US
    W::WS
    obstacles::OB
    f_eval::FE
    f_backward_eval::FB
    Uformat::UF
    Wformat::WF
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
