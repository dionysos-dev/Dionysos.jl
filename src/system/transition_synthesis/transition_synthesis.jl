# Ellipsoid-to-ellipsoid transition synthesis for affine systems: given
# x⁺ = Ax + Bu + c (+ Dw), decide whether an affine controller
# u(x) = K(x − c₁) + ℓ drives a source ellipsoid into a target ellipsoid under
# input constraints and polytopic noise, while minimizing an upper bound J on
# the worst-case transition cost ‖Λ·[x; u; 1]‖². Each entry point assembles one
# SDP from shared S-procedure blocks; see Corollary 1 of
# https://arxiv.org/pdf/2204.00315.pdf.
#
# Entry points:
#   `solve_transition`                — both ellipsoids fixed, controller free;
#   `solve_transition_backward`       — target fixed, source shape synthesized;
#   `solve_transition_backward_2step` — two-step source synthesis through a
#                                       fixed second-step controller;
#   `solve_transition_forward`        — source fixed, target synthesized.
#
# Layout (one concern per file):
#   types.jl        — TransitionResult, AffineSys, set → LMI-data normalization;
#   blocks.jl       — the S-procedure PSD block builders (pure matrix functions);
#   kernel_tools.jl — the shared kernel skeleton: pose-and-record constraints
#                     (one argument tuple poses the JuMP constraint AND rebuilds
#                     the numeric block for a-posteriori validation), remainder
#                     models, status predicate, size objectives;
#   fixed.jl / backward.jl / two_step.jl / forward.jl — one kernel each.

import LazySets
using JuMP

include("types.jl")
include("blocks.jl")
include("kernel_tools.jl")
include("fixed.jl")
include("backward.jl")
include("two_step.jl")
include("forward.jl")
