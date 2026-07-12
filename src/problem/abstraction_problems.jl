# ----------------------------
# Abstraction problems
# ----------------------------

"""
    AlternatingSimulationProblem{S, X} <: AbstractionProblem

A problem type used to construct a **sound abstraction** of a dynamical system.

- `S`: The system to abstract (continuous or discrete-time).
- `X`: The state-space region of interest to abstract.

This problem encodes no control objective. It is intended for generating symbolic models that can later be reused by other solvers.
"""
struct AlternatingSimulationProblem{S, X} <: AbstractionProblem
    system::S
    state_set::X
end

"""
    BisimulationQuotientProblem{S, X, R} <: AbstractionProblem

A problem type used to construct a finite bisimulation (exact equivalence abstraction) quotient induced.

# Fields
- `system`: the switched system to abstract.
- `state_set`: the state-space region of interest `X`.
- `observation_regions`: the regions of interest used to define the observation map.

This problem encodes no control objective. It is intended for generating symbolic models that can later be reused by other solvers.
"""
struct BisimulationQuotientProblem{S, X, R} <: AbstractionProblem
    system::S
    state_set::X
    observation_regions::R
end
