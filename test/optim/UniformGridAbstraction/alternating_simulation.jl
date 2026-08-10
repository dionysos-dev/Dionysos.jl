module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI
import Suppressor
import LazySets

include("../../../problems/Integrator/integrator.jl")

# Directly exercises the abstraction-only optimizer `OptimizerAlternatingSimulationProblem`
# across its approximation modes and mapping options. Only the GROWTH continuous path is hit
# by the higher-level solver tests; the other modes and the implicit/periodic mapping paths
# are covered here.
const UGA = AB.UniformGridAbstraction

_attr(name) = MOI.RawOptimizerAttribute(name)

function _set!(opt, pairs)
    for (name, value) in pairs
        MOI.set(opt, _attr(name), value)
    end
    return opt
end

# Build the abstraction and return the resulting symbolic model. `@suppress` hides the
# per-build "noise not accounted for" warning.
function _build(opt)
    return Suppressor.@suppress begin
        MOI.optimize!(opt)
        MOI.get(opt, _attr("abstract_system"))
    end
end

function _base_optimizer(asp; approx_mode)
    opt = UGA.OptimizerAlternatingSimulationProblem()
    _set!(
        opt,
        [
            "alternating_simulation_problem" => asp,
            "input_grid" => MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5)),
            "approx_mode" => approx_mode,
            "print_level" => 0,
        ],
    )
    return opt
end

const _STATE_GRID = MP.GridFree(SVector(-2.0, -2.0), SVector(0.5, 0.5))

@testset "AlternatingSimulation: continuous approx modes" begin
    csys = Integrator.system()
    asp = PR.AlternatingSimulationProblem(csys, csys.X)

    # CENTER_SIMULATION: no extra data required.
    opt = _base_optimizer(asp; approx_mode = UGA.CENTER_SIMULATION)
    _set!(opt, ["state_grid" => _STATE_GRID, "time_step" => 0.3])
    @test SY.get_n_state(_build(opt)) > 0

    # RANDOM_SIMULATION: sample a few points per cell.
    opt = _base_optimizer(asp; approx_mode = UGA.RANDOM_SIMULATION)
    _set!(opt, ["state_grid" => _STATE_GRID, "time_step" => 0.3, "n_samples" => 3])
    @test SY.get_n_state(_build(opt)) > 0

    # LINEARIZED: ẋ = u has a zero Jacobian and Hessian.
    opt = _base_optimizer(asp; approx_mode = UGA.LINEARIZED)
    _set!(
        opt,
        [
            "state_grid" => _STATE_GRID,
            "time_step" => 0.3,
            "DF_sys" => (x, u) -> (@SMatrix zeros(2, 2)),
            "bound_DF" => u -> 0.0,
            "bound_DDF" => u -> 0.0,
        ],
    )
    @test SY.get_n_state(_build(opt)) > 0

    # GROWTH via a directly-supplied growth-bound map. A continuous growth-bound map takes
    # `(radius, u, tstep)`; for ẋ = u the radius is unchanged.
    opt = _base_optimizer(asp; approx_mode = UGA.GROWTH)
    _set!(
        opt,
        [
            "state_grid" => _STATE_GRID,
            "time_step" => 0.3,
            "growthbound_map" => (r, u, tstep) -> r,
        ],
    )
    @test SY.get_n_state(_build(opt)) > 0
end

@testset "AlternatingSimulation: discrete approx modes" begin
    f(x, u) = x .+ 0.3 .* u
    dsys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        f,
        2,
        2,
        LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0)),
        LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0)),
    )
    asp = PR.AlternatingSimulationProblem(dsys, dsys.X)

    # A discrete-time system skips the continuous discretization path.
    opt = _base_optimizer(asp; approx_mode = UGA.CENTER_SIMULATION)
    _set!(opt, ["state_grid" => _STATE_GRID])
    @test SY.get_n_state(_build(opt)) > 0

    opt = _base_optimizer(asp; approx_mode = UGA.GROWTH)
    _set!(opt, ["state_grid" => _STATE_GRID, "growthbound_map" => (r, u) -> r])
    @test SY.get_n_state(_build(opt)) > 0

    opt = _base_optimizer(asp; approx_mode = UGA.RANDOM_SIMULATION)
    _set!(opt, ["state_grid" => _STATE_GRID, "n_samples" => 3])
    @test SY.get_n_state(_build(opt)) > 0
end

@testset "AlternatingSimulation: implicit mapping and state set" begin
    csys = Integrator.system()
    asp = PR.AlternatingSimulationProblem(csys, csys.X)

    opt = _base_optimizer(asp; approx_mode = UGA.CENTER_SIMULATION)
    _set!(
        opt,
        [
            "state_grid" => _STATE_GRID,
            "time_step" => 0.3,
            "use_implicit_mapping" => true,
            "use_implicit_stateset" => true,
            # must enclose the abstraction region (the system domain [-2, 2]²)
            "mapping_region" => LazySets.Hyperrectangle(;
                low = SVector(-2.5, -2.5),
                high = SVector(2.5, 2.5),
            ),
        ],
    )
    @test SY.get_n_state(_build(opt)) > 0
end

@testset "AlternatingSimulation: periodic mapping built from h" begin
    csys = Integrator.system()   # domain [-2, 2]² — width 4 in each dim
    asp = PR.AlternatingSimulationProblem(csys, csys.X)

    # No explicit state_grid: the grid is built from `h`, periodic in dim 1 with period 4.
    opt = _base_optimizer(asp; approx_mode = UGA.CENTER_SIMULATION)
    _set!(
        opt,
        [
            "h" => SVector(0.5, 0.5),
            "time_step" => 0.3,
            "use_periodic_mapping" => true,
            "periodic_dims" => SVector(1),
            "periodic_periods" => SVector(4.0),
            "periodic_start" => SVector(-2.0),
        ],
    )
    @test SY.get_n_state(_build(opt)) > 0
end

end # module TestMain
