module Abstraction

include("UniformGridAbstraction/uniform_grid_abstraction.jl")
include("UniformEllipsoidAbstraction/uniform_ellipsoid_abstraction.jl")
include("HybridSystemAbstraction/hybrid_system_abstraction.jl")
include("lazy_ellipsoids_abstraction.jl")

include("PathCompleteBisimulation/quotient_bisimulation_problem.jl")
end
