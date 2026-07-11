# Local affine approximation of nonlinear dynamics: an affine system valid around a
# linearization point together with Lipschitz bounds on the linearization error.
# The symbolic implementation lives in the Symbolics package extension; the stubs
# below error until it is loaded.

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

# ------------------------------------------------------------
# Affine-approximation providers (package extension hooks)
# ------------------------------------------------------------
export AbstractAffineApproximationProvider,
    SymbolicAffineApproximationProvider, AffineApproximation, build_affine_approximation

abstract type AbstractAffineApproximationProvider end

struct SymbolicAffineApproximationProvider{F, X, U, W, DW, UF, WF} <:
       AbstractAffineApproximationProvider
    fsymbolic::F
    x::X
    u::U
    w::W
    ΔW::DW
    Uformat::UF
    Wformat::WF
end

struct AffineApproximation{SYS, LIP, UF, WF, S}
    system::SYS
    lipschitz::LIP
    Uformat::UF
    Wformat::WF
    summary::S
end

function build_affine_approximation(
    provider::AbstractAffineApproximationProvider,
    k::Int,
    xk,
    xnext,
    uk,
    δx,
    δu,
)
    return error("build_affine_approximation not implemented for $(typeof(provider))")
end
