# Local affine approximation of nonlinear dynamics: an affine system valid around a
# linearization point together with Lipschitz bounds on the linearization error.
# One entry point for every consumer (lazy-ellipsoids abstraction, ellipsoidal
# backward certifier):
#
#     build_affine_approximation(provider, x̄, ū, w̄ = nothing; δx, δu) -> AffineApproximation
#
# Providers: `SymbolicAffineApproximationProvider` (implemented in the Symbolics
# package extension) and `AnalyticAffineApproximationProvider` (user-supplied
# derivative information, no extension needed).

export AbstractAffineApproximationProvider,
    SymbolicAffineApproximationProvider,
    AnalyticAffineApproximationProvider,
    AffineApproximation,
    build_affine_approximation

"""
    AbstractAffineApproximationProvider

Supertype of the strategies producing local affine approximations of the
dynamics; see [`build_affine_approximation`](@ref).
"""
abstract type AbstractAffineApproximationProvider end

"""
    AffineApproximation

Result of [`build_affine_approximation`](@ref): the local affine `system`
(`MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem`), the
`lipschitz` bounds on the linearization error over the local domain, the
`Uformat`/`Wformat` matrices consumed by the transition-synthesis LMIs, and a
provider-specific `summary` NamedTuple.
"""
struct AffineApproximation{SYS, LIP, UF, WF, S}
    system::SYS
    lipschitz::LIP
    Uformat::UF
    Wformat::WF
    summary::S
end

"""
    build_affine_approximation(provider, x̄, ū, w̄ = nothing; δx, δu) -> AffineApproximation

Affine approximation of the dynamics around the linearization point `(x̄, ū, w̄)`,
valid on the box of radii `δx` / `δu` around it (`w̄ = nothing` linearizes at zero
noise). Implemented per provider.
"""
function build_affine_approximation(
    provider::AbstractAffineApproximationProvider,
    x̄,
    ū,
    w̄ = nothing;
    δx,
    δu,
)
    return error("build_affine_approximation not implemented for $(typeof(provider))")
end

# ------------------------------------------------------------
# Symbolic provider (implementation in the Symbolics extension)
# ------------------------------------------------------------

"""
    SymbolicAffineApproximationProvider

Affine approximation from symbolic dynamics `fsymbolic(x, u, w)`: Jacobians are
computed symbolically and the Lipschitz constants are bounded over the local
domain with interval arithmetic. `ΔW` is the noise-box radius vector; `Uformat` /
`Wformat` are the LMI encodings of the input and noise sets
([`format_input_set`](@ref) / [`format_noise_set`](@ref)). Requires Symbolics.jl
(`using Symbolics`).
"""
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

# ------------------------------------------------------------
# Analytic provider (no extension needed)
# ------------------------------------------------------------

"""
    AnalyticAffineApproximationProvider(; A, B, f, lipschitz, nw = 1, E = nothing,
                                        ΔW = zeros(nw), Uformat, Wformat)

Affine approximation from user-supplied derivative information — the same data
UGA's `LINEARIZED` mode asks for, usable without Symbolics:

- `A(x̄, ū, w̄)`, `B(x̄, ū, w̄)`: state and input Jacobians at the linearization point;
- `E(x̄, ū, w̄)`: noise matrix (defaults to a zero `nx × nw` matrix);
- `f(x̄, ū, w̄)`: the dynamics value, giving the affine offset `c = f − A·x̄ − B·ū − E·w̄`;
- `lipschitz`: error bound over the local box — either a constant vector or a
  callable `(x̄, ū, w̄, δx, δu) -> Vector`;
- `Uformat` / `Wformat`: as in [`SymbolicAffineApproximationProvider`](@ref).
"""
struct AnalyticAffineApproximationProvider{FA, FB, FE, FF, L, DW, UF, WF} <:
       AbstractAffineApproximationProvider
    A::FA
    B::FB
    E::FE
    f::FF
    lipschitz::L
    nw::Int
    ΔW::DW
    Uformat::UF
    Wformat::WF
end

function AnalyticAffineApproximationProvider(;
    A,
    B,
    f,
    lipschitz,
    nw::Int = 1,
    E = nothing,
    ΔW = zeros(nw),
    Uformat,
    Wformat,
)
    Efun = E === nothing ? ((x̄, ū, w̄) -> zeros(length(x̄), nw)) : E
    return AnalyticAffineApproximationProvider(
        A,
        B,
        Efun,
        f,
        lipschitz,
        nw,
        ΔW,
        Uformat,
        Wformat,
    )
end

_lipschitz_bound(L::AbstractVector, x̄, ū, w̄, δx, δu) = L
_lipschitz_bound(L, x̄, ū, w̄, δx, δu) = L(x̄, ū, w̄, δx, δu)

function build_affine_approximation(
    provider::AnalyticAffineApproximationProvider,
    x̄,
    ū,
    w̄ = nothing;
    δx,
    δu,
)
    w̄ === nothing && (w̄ = zeros(provider.nw))

    A = provider.A(x̄, ū, w̄)
    B = provider.B(x̄, ū, w̄)
    E = provider.E(x̄, ū, w̄)
    c = provider.f(x̄, ū, w̄) - A * collect(x̄) - B * collect(ū) - E * collect(w̄)

    Xbar = LazySets.Hyperrectangle(; low = collect(x̄) .- δx, high = collect(x̄) .+ δx)
    Ubar = LazySets.Hyperrectangle(; low = collect(ū) .- δu, high = collect(ū) .+ δu)
    Wbar = LazySets.Hyperrectangle(;
        low = collect(w̄) .- provider.ΔW,
        high = collect(w̄) .+ provider.ΔW,
    )

    system = MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, Xbar, Ubar, Wbar)
    L = _lipschitz_bound(provider.lipschitz, x̄, ū, w̄, δx, δu)

    return AffineApproximation(
        system,
        L,
        provider.Uformat,
        provider.Wformat,
        (; δx = copy(δx), δu = copy(δu)),
    )
end

# ------------------------------------------------------------
# Symbolic problem bundle (consumed by the lazy-ellipsoids solver)
# ------------------------------------------------------------

"""
    SymbolicSystem

System-definition bundle for solvers that re-linearize the dynamics on the fly
(lazy-ellipsoids abstraction): symbolic dynamics `(fsymbolic, x, u, w)`,
dimensions, domains `X`/`U`/`W`, local linearization radii `ΔX`/`ΔU`/`ΔW`,
forward/backward evaluators, and the LMI formats of `U` and `W`. Use
[`get_affine_provider`](@ref) to obtain the corresponding
[`SymbolicAffineApproximationProvider`](@ref).
"""
struct SymbolicSystem{FS, XT, UTT, WT, DX, DU, DW, XS, US, WS, FE, FB, UF, WF}
    fsymbolic::FS
    nx::Int
    nu::Int
    nw::Int
    x::XT
    u::UTT
    w::WT
    ΔX::DX
    ΔU::DU
    ΔW::DW
    X::XS
    U::US
    W::WS
    f_eval::FE
    f_backward_eval::FB
    Uformat::UF
    Wformat::WF
end

"""
    get_affine_provider(sys::SymbolicSystem) -> SymbolicAffineApproximationProvider

The affine-approximation provider matching the bundle's symbolic dynamics.
"""
function get_affine_provider(sys::SymbolicSystem)
    return SymbolicAffineApproximationProvider(
        sys.fsymbolic,
        sys.x,
        sys.u,
        sys.w,
        sys.ΔW,
        sys.Uformat,
        sys.Wformat,
    )
end
