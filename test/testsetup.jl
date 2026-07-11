# Shared test harness — included by every test file so they share one set of imports,
# module aliases, tiny fixtures and assertion helpers instead of re-declaring them.
#
# Usage (works both under runtests.jl and when a file is run standalone):
#
#     using Dionysos
#     include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))
#
# The path is resolved from `pathof(Dionysos)` so it is independent of the caller's depth.
# Everything below is defined in the *includer's* module scope, keeping each test file a
# self-contained, standalone-runnable module.

using Test
using StaticArrays
using MathematicalSystems

# Standard module aliases (match src/Dionysos.jl and the house style).
const DI = Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems
const AB = OP.Abstraction

# ----------------------------------------------------------------------------------------
# Tiny canonical systems — deterministic and small so abstractions are exact and fast.
# ----------------------------------------------------------------------------------------

"""
    single_integrator(; n = 1, xbound = 2.0, ubound = 1.0)

`ẋ = u` on `[-xbound, xbound]^n` with inputs in `[-ubound, ubound]^n`.
The simplest controllable test plant; drive to a target and stay with `u = 0`.
"""
function single_integrator(; n::Int = 1, xbound = 2.0, ubound = 1.0)
    f(x, u) = u
    X = UT.box(SVector{n}(fill(-xbound, n)...), SVector{n}(fill(xbound, n)...))
    U = UT.box(SVector{n}(fill(-ubound, n)...), SVector{n}(fill(ubound, n)...))
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(f, n, n, X, U)
end

"""
    double_integrator(; xbound = 2.0, vbound = 2.0, ubound = 1.0)

`ẋ = v, v̇ = u` (2D state, 1D input): a minimal second-order plant for reach/safety tests.
"""
function double_integrator(; xbound = 2.0, vbound = 2.0, ubound = 1.0)
    f(x, u) = SVector(x[2], u[1])
    X = UT.box(SVector(-xbound, -vbound), SVector(xbound, vbound))
    U = UT.box(SVector(-ubound), SVector(ubound))
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(f, 2, 1, X, U)
end

# ----------------------------------------------------------------------------------------
# Assertion helpers
# ----------------------------------------------------------------------------------------

"""
    controller_admissible(ctrl, x) -> Bool

`true` when `ctrl` is defined at state `x` and yields a control there. Use inside `@test`.
"""
function controller_admissible(ctrl, x)
    ST.is_defined(ctrl, nothing, x) || return false
    return ST.output_control(ctrl, nothing, x) !== nothing
end
