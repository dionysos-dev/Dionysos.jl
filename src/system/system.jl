module System

import MathematicalSystems as MS
import HybridSystems

import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA

import JuMP: MOI
import RecipesBase: @recipe, @series

using ..Utils
UT = Utils

# Continuous systems: MathematicalSystems extension
include("continuous_systems/continuous_systems.jl")

# Time discretization
include("time_discretization.jl")

# Linearization
include("linearization.jl")

# Kernel Approximation
include("kernel_approximation/kernel_approximation.jl")

# Controllers
include("controller.jl")
include("pid_controller.jl")

# Trajectories
include("trajectories.jl")

# ------------------------------------------------------------
# Package extension hooks
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

end
