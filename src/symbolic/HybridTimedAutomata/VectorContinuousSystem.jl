
@warn "Ce type de système mathématique est nécessaire pour réaliser le produit cartésien entre le système dynamique et le temps. À terme, ce code devrait être intégré à MathematicalSystems.jl et ne devrait plus se trouver dans le code de Dionysos, ou du moins pas dans la partie symbolic."

"""
    VectorContinuousSystem(systems::Vector{<:AbstractContinuousSystem})

Container for a vector of mathematical systems, allowing joint manipulation as a single system.
"""
struct VectorContinuousSystem <: AbstractContinuousSystem
    systems::Vector{AbstractContinuousSystem}
end

#pas défini pour les ConstrainedBlackBoxControlContinuousSystem et pas particulièrement utile pour le moment 
stateset(sys::VectorContinuousSystem) = tuple((MathematicalSystems.stateset(s) for s in sys.systems)...)
inputset(sys::VectorContinuousSystem) = tuple((MathematicalSystems.inputset(s) for s in sys.systems)...)
statedim(sys::VectorContinuousSystem) = sum(MathematicalSystems.statedim(s) for s in sys.systems)
inputdim(sys::VectorContinuousSystem) = sum(MathematicalSystems.inputdim(s) for s in sys.systems)
