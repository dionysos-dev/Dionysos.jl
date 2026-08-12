module System

import MathematicalSystems as MS
import HybridSystems

import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA

import JuMP: MOI
import RecipesBase: @recipe, @series
import LazySets

using ..Utils
UT = Utils

# Continuous systems: MathematicalSystems extension
include("vector_continuous_system.jl")

# Hybrid transitions: the guard/reset pair stored on a HybridSystems transition
include("reset_map.jl")

# Time discretization
include("time_discretization.jl")

# Local affine approximation of nonlinear dynamics (linearization + error bounds)
include("affine_approximation.jl")

# Reachable-set over/under-approximations of the dynamics
include("approximation/approximation.jl")

# Controllers
include("controllers/controller.jl")
include("controllers/pid_controller.jl")

# Trajectories & closed-loop simulation
include("trajectories/trajectory.jl")
include("trajectories/closed_loop.jl")

# Normalized certification frame z = x ./ t (set/trajectory/controller half;
# the dynamics half is `normalized_symbolic_provider` in the Symbolics ext)
include("normalized_frame.jl")

# Synthesis of controlled ellipsoid-to-ellipsoid transitions (S-procedure LMIs)
include("transition_synthesis/transition_synthesis.jl")

# Lyapunov state feedback (stabilizability LMI)
include("stabilizing_feedback.jl")

end
