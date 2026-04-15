module System

import MathematicalSystems as MS
import HybridSystems

import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA

import JuMP: MOI
import RecipesBase: @recipe, @series

using ..Utils
UT = Utils

# Discrete systems: HybridSystems extension
include("discrete_systems/discrete_systems.jl")

# Continuous systems: MathematicalSystems extension
include("continuous_systems/continuous_systems.jl")

# Time discretization
include("time_discretization.jl")

# Linearization
include("linearization.jl")

# Kernel Approximation
include("kernel_approximation/kernel_approximation.jl")

# Basic Controllers
include("pid_controller.jl")

# Trajectories
include("trajectories.jl")

end
