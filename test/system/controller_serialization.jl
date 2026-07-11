# The user-level guarantee behind the "controllers are plain data" convention:
# every synthesized controller can be saved to JLD2 and reloaded in a fresh
# session, and the reloaded controller produces the same controls.
module TestControllerSerialization

using Test
using StaticArrays
using JLD2
import MathematicalSystems as MS
import MathOptInterface as MOI

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

function roundtrip(obj)
    path = joinpath(mktempdir(), "controller.jld2")
    jldsave(path; obj = obj)
    return jldopen(path, "r") do f
        return f["obj"]
    end
end

@testset "Table-backed controllers round-trip" begin
    tab = ST.ControlTable(3)
    ST.add_control!(tab, 1, 4)
    ST.set_control!(tab, 2, 7)
    ctrl = ST.DiscreteStaticController(Set([1, 2]), tab, false)

    ctrl2 = roundtrip(ctrl)
    for q in 1:3
        @test ST.output_control(ctrl2, nothing, q) == ST.output_control(ctrl, nothing, q)
    end

    affine = ST.AffineController(MS.AffineMap([1.0 0.0; 0.0 2.0], [0.5, -0.5]))
    affine2 = roundtrip(affine)
    x = [1.0, -1.0]
    @test ST.output_control(affine2, nothing, x) == ST.output_control(affine, nothing, x)

    memtab = ST.ControlTable(2)
    ST.add_control!(memtab, 1, 9)
    ST.add_control!(memtab, 2, 8)
    mem = ST.AutomatonMemoryController(
        1,
        [2, 1],
        1,
        Dict((1, 1) => 1, (1, 2) => 2, (2, 1) => 2, (2, 2) => 2),
        Dict((1, 1) => 1, (2, 2) => 2),
        ST.DiscreteStaticController(Set([1, 2]), memtab, false),
    )
    mem2 = roundtrip(mem)
    @test ST.initial_state(mem2) == ST.initial_state(mem)
    for qa in 1:2, qs in 1:2
        @test ST.output_control(mem2, qa, qs) == ST.output_control(mem, qa, qs)
        @test ST.update_state(mem2, qa, qs) == ST.update_state(mem, qa, qs)
    end
end

@testset "Synthesized UGA controller round-trip" begin
    # Small reachability problem end to end, then save/load the concrete controller.
    F_sys(x, u) = SVector(u[1], -0.5 * x[1])
    _X_ = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    _U_ = UT.box(SVector(-1.0), SVector(1.0))
    concrete_system = MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 2, 1, _X_, _U_)

    target = UT.box(SVector(0.5, -0.5), SVector(1.5, 0.5))
    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        UT.box(SVector(-1.5, -0.5), SVector(-0.5, 0.5)),
        target,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0, 0.0), SVector(0.25, 0.25)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.25)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.4)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        u -> SMatrix{2, 2}(0.0, 0.5, 0.0, 0.0),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    ctrl = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    @test ctrl !== nothing

    ctrl2 = roundtrip(ctrl)

    npts = 0
    for x1 in -1.9:0.4:1.9, x2 in -1.9:0.4:1.9
        x = SVector(x1, x2)
        u = ST.output_control(ctrl, nothing, x)
        u2 = ST.output_control(ctrl2, nothing, x)
        @test u == u2
        u === nothing || (npts += 1)
    end
    @test npts > 0    # the comparison exercised actual controls, not only nothings
end

end # module
