# """
#     VectorContinuousSystem{T<:AbstractContinuousSystem}

# Container for a vector of mathematical systems, allowing joint manipulation as a single system.
# This type performs the Cartesian product between dynamical systems and time components.

# # Note
# Eventually, this code should be integrated into MathematicalSystems.jl and should no longer 
# be present in the Dionysos codebase.

# # Fields
# - `systems::Vector{T}`: Vector of continuous systems to be composed

# # Example
# ```julia
# dyn_sys = ConstrainedBlackBoxControlContinuousSystem(...)
# time_sys = ConstrainedLinearContinuousSystem(...)
# composed = VectorContinuousSystem([dyn_sys, time_sys])
# ```
# """
struct VectorContinuousSystem{T <: AbstractContinuousSystem} <: AbstractContinuousSystem
    systems::Vector{T}
end

# Constructor for mixed types
VectorContinuousSystem(systems::Vector{<:AbstractContinuousSystem}) =
    VectorContinuousSystem{AbstractContinuousSystem}(systems)

MathematicalSystems.stateset(sys::VectorContinuousSystem) =
    Tuple(MathematicalSystems.stateset(s) for s in sys.systems)

MathematicalSystems.inputset(sys::VectorContinuousSystem) =
    Tuple(MathematicalSystems.inputset(s) for s in sys.systems)

MathematicalSystems.statedim(sys::VectorContinuousSystem)::Int =
    sum(MathematicalSystems.statedim(s) for s in sys.systems; init = 0)

MathematicalSystems.inputdim(sys::VectorContinuousSystem)::Int =
    sum(MathematicalSystems.inputdim(s) for s in sys.systems; init = 0)
