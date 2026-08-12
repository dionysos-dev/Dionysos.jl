module EllipsoidalTrajectoryCertifier

# Ellipsoidal funnel certification of a nominal trajectory (plan.md §4–§5): one
# transition-synthesis SDP per step chains ellipsoids (E_1, κ_1, …, E_{K+1}) with the
# guarantee x ∈ E_k, u = κ_k(x) ⟹ x⁺ ∈ E_{k+1}. The module is split by concern:
#
#   options.jl        — ChainOptions + AdaptiveLinearizationBoxOptions
#   diagnostics.jl    — step summaries, radii/volume helpers (out of the hot path)
#   context.jl        — ChainContext (problem + trajectory + provider), state scaling
#   steps.jl          — StepRecord and the fixed-box backward step
#   adaptive_boxes.jl — the adaptive linearization-box search (backward only)
#   two_step.jl       — the two-step rescue (skip a needle intermediate funnel)
#   gates.jl          — soundness gates: box consistency, collapse, state domain,
#                       terminal containment, initial coverage
#   chain.jl          — run_chain! + CertificationResult
#   certifier.jl      — BackwardCertifier (the AbstractTrajectoryCertifier front)

import Dionysos
import LinearAlgebra as LA
import MathematicalSystems as MS
import LazySets

const DI = Dionysos
const PR = DI.Problem
const UT = DI.Utils
const ST = DI.System

import ..AbstractTrajectoryCertifier
import ..set_problem!
import ..set_trajectory!
import ..certify!
import ..get_controller
import ..get_success
import ..get_solve_time
# Generator half of the interface, used by the re-planning loop.
import ..set_seed_trajectory!
import ..generate!
import ..get_trajectory

export BackwardCertifier,
    ForwardCertifier,
    AdaptiveLinearizationBoxOptions,
    ChainOptions,
    ForwardOptions,
    StepRecord,
    CertificationResult,
    get_result,
    build_context,
    backward_step!,
    forward_step!,
    run_chain!,
    bidirectional_certify!,
    prefix_replan_certify!

include("options.jl")
include("diagnostics.jl")
include("context.jl")
include("steps.jl")
include("adaptive_boxes.jl")
include("two_step.jl")
include("gates.jl")
include("chain.jl")
include("certifier.jl")
include("forward.jl")
include("bidirectional.jl")
include("replan.jl")

end # module
