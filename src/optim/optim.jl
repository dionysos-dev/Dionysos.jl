module Optim

include("optimizer_common.jl")
include("lifted_control_optimizer.jl")

include("discrete_systems/discrete_systems.jl")

module Abstraction
    include("continuous_systems/continuous_systems.jl")
    include("hybrid_systems/hybrid_systems.jl")
    include("trajectory_generators/trajectory_generators.jl")
    include("trajectory_certifiers/trajectory_certifiers.jl")
end

include("bemporad_morari.jl")
include("branch_and_bound.jl")
# Q-function / dual-dynamic-programming lower-bound algorithms consumed by BranchAndBound
# (not reinforcement Q-learning).
include("q_function.jl")

end
