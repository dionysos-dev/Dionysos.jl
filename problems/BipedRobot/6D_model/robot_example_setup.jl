# ==============================================================================
#  Setup module for robot_example — problem definition, grids, abstract system.
#  No addprocs, no @everywhere. Used by partition_runner.jl, assemble_partitions.jl,
#  and (optionally) refactored robot_example.jl.
# ==============================================================================

module RobotExampleSetup

using MathematicalSystems
using StaticArrays
using LinearAlgebra
using JuMP

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(@__DIR__, "robot_problem.jl"))
include(joinpath(@__DIR__, "utils.jl"))

const ROBOT_URDF = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")
const TSTEP = 0.1

"""
    setup(; partition_strategy, state_filter, state_input_filter)

Build the concrete system, grids, empty abstract system (no transitions),
and system approximation.  Returns a NamedTuple with everything a partition
runner or assembler needs.

The abstract system construction replicates the internal steps of
`OptimizerEmptyProblem` so that the state/input numbering is identical
to a standard `MOI.optimize!` run.
"""
function setup(;
    partition_strategy::Symbol = :roundrobin,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    # ---- Concrete system --------------------------------------------------
    concrete_problem = RobotProblem.problem(; robot_urdf = ROBOT_URDF, tstep = TSTEP)
    concrete_system  = concrete_problem.system

    n_state = MathematicalSystems.statedim(concrete_system)
    n_input = MathematicalSystems.inputdim(concrete_system)

    # ---- Grids (same as robot_example.jl) ---------------------------------
    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx = SVector{n_state, Float64}([fill(2π / 180, 3)..., fill(0.15, 3)...]) * 2.5
    state_grid = MP.GridFree(x0, hx)

    u0 = SVector{n_input, Float64}(zeros(n_input))
    hu = SVector{n_input, Float64}(fill(1.0, n_input))
    input_grid = MP.GridFree(u0, hu)

    # ---- Replicate optimizer setup (empty_problem.jl MOI.optimize!) -------
    # Abstraction region  (same as build_abstraction_region)
    abstraction_region = concrete_system.X      # HyperRectangle

    # State mapping  (same as build_state_mapping with defaults)
    N = MP.get_dim(state_grid)
    XMapping = MP.ExplicitGridMapping{N, Float64}(state_grid)
    MP.add_set!(XMapping, abstraction_region, MP.INNER)

    # Input mapping  (same as build_input_mapping)
    M = MP.get_dim(input_grid)
    UMapping = MP.ExplicitGridMapping{M, Float64}(
        input_grid, concrete_system.U, MP.CENTER,
    )

    # State set  (same as build_state_set with defaults)
    Xset = MP.ExplicitIdSet{N}()
    MP.add_set!(Xset, XMapping, abstraction_region, MP.INNER)

    # Retained set  (same as build_allowed_state_set)
    Rset = copy(Xset)

    # Abstract system (empty — no transitions yet)
    abstract_system = SY.SymbolicModelList(
        XMapping, UMapping;
        Xset = Xset,
        Rset = Rset,
        # Uset defaults to MappingSet{M}() inside the constructor
    )

    # System approximation  (CENTER_SIMULATION, efficient = true)
    system_approximation = ST.DiscreteTimeCenteredSimulation(concrete_system)

    return (
        concrete_problem      = concrete_problem,
        concrete_system       = concrete_system,
        abstract_system       = abstract_system,
        system_approximation  = system_approximation,
        state_grid            = state_grid,
        input_grid            = input_grid,
        n_state               = n_state,
        n_input               = n_input,
        partition_strategy    = partition_strategy,
        state_filter          = state_filter,
        state_input_filter    = state_input_filter,
    )
end

end # module RobotExampleSetup
