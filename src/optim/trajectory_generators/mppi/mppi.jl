module MPPITrajectoryGenerator

# MPPI as a one-shot trajectory optimizer (plan.md §3): softmin/:cem updates over
# threaded, penalty-scored rollouts with ESS-adaptive temperature, structured noise,
# and an optional importance-sampling correction.

import MathematicalSystems as MS
import Dionysos
using StaticArrays

import ..AbstractTrajectoryGenerator
import ..set_problem!
import ..generate!
import ..get_trajectory
import ..get_success
import ..get_solve_time
import ..AbstractCostTerm
import ..CompositeCost
import ..rollout_cost
import ..rollout_trajectory

const DI = Dionysos
const ST = DI.System
const PR = DI.Problem

export GaussianMPPINoise

include("noise.jl")
include("updates.jl")
include("generator.jl")

end # module
