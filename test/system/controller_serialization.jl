# The user-level guarantee behind the "controllers are plain data" convention:
# every synthesized controller can be saved to JLD2 and reloaded in a fresh
# session, and the reloaded controller produces the same controls.
module TestControllerSerialization

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using JLD2
import MathematicalSystems as MS
import MathOptInterface as MOI
import LazySets

function roundtrip(obj)
    path = joinpath(mktempdir(), "controller.jld2")
    jldsave(path; obj = obj)
    return jldopen(path, "r") do f
        return f["obj"]
    end
end

"""
    fresh_session_answers(ctrl, states) -> (ok, err, answers)

Save `ctrl`, then load it in a **separate Julia process** and collect
`(is_defined, output_control)` on each of `states`.

Reloading in the same process is not the guarantee we care about and cannot catch
the failure that matters: a field whose type no longer exists on load — a closure
is the usual culprit — still resolves in-process, because the type is defined in
the running session. Deployment loads the file in a fresh session, where JLD2
substitutes a non-callable placeholder and the controller breaks. So the check has
to cross a process boundary.
"""
function fresh_session_answers(ctrl, states)
    dir = mktempdir()
    ctrl_path = joinpath(dir, "controller.jld2")
    states_path = joinpath(dir, "states.jld2")
    out_path = joinpath(dir, "answers.jld2")
    jldsave(ctrl_path; controller = ctrl)
    jldsave(states_path; states = states)

    # `global` is load-bearing: a `try` body is a soft scope, so a plain assignment
    # here would bind a local and leave the outer values untouched — the child would
    # then report success with no answers.
    code = """
    import Dionysos, JLD2
    const ST = Dionysos.System
    ok, err, answers = true, "", Any[]
    try
        ctrl = JLD2.load(raw"$(ctrl_path)", "controller")
        states = JLD2.load(raw"$(states_path)", "states")
        mem = ST.initial_state(ctrl)
        global answers = Any[(ST.is_defined(ctrl, mem, x),
                              ST.output_control(ctrl, mem, x)) for x in states]
    catch e
        global ok = false
        global err = sprint(showerror, e)
    end
    JLD2.jldsave(raw"$(out_path)"; ok = ok, err = err, answers = answers)
    """

    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $code`
    success(pipeline(cmd; stdout = devnull, stderr = devnull)) ||
        return (false, "the child process itself failed", Any[])

    return jldopen(out_path, "r") do f
        return (f["ok"], f["err"], f["answers"])
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
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
    concrete_system = MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 2, 1, _X_, _U_)

    target = LazySets.Hyperrectangle(; low = SVector(0.5, -0.5), high = SVector(1.5, 0.5))
    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        LazySets.Hyperrectangle(; low = SVector(-1.5, -0.5), high = SVector(-0.5, 0.5)),
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
        # `is_defined` is checked as well as `output_control`: it is the method that
        # reads the controller's *domain*, so a field that fails to deserialize can
        # leave a controller still emitting commands while silently having lost the
        # "is this state controlled at all?" guarantee.
        @test ST.is_defined(ctrl2, nothing, x) == ST.is_defined(ctrl, nothing, x)
        u = ST.output_control(ctrl, nothing, x)
        u2 = ST.output_control(ctrl2, nothing, x)
        @test u == u2
        u === nothing || (npts += 1)
    end
    @test npts > 0    # the comparison exercised actual controls, not only nothings
end

@testset "Slew-rate limited UGA controller round-trip" begin
    # `BoundedInputVariation` is the controller shape that gets deployed
    # (`control_server`), and the constraint is the one place the pipeline stores a
    # *function* rather than plain data — so it is where the "controllers are plain
    # data" guarantee is most likely to break.
    #
    # It also requires a *deterministic* abstraction, so this follows the same recipe
    # as the biped example: `ẋ = u` on a lattice where `tstep * du == dx`, abstracted
    # with `CENTER_SIMULATION`, which makes the abstraction exact (a bisimulation).
    dx, du, tstep = 0.25, 0.5, 0.5    # tstep * du == dx

    F_sys(x, u) = SVector(u[1])
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0), high = SVector(2.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
    concrete_system = MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 1, 1, _X_, _U_)

    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        LazySets.Hyperrectangle(; low = SVector(-1.25), high = SVector(-0.75)),
        LazySets.Hyperrectangle(; low = SVector(0.75), high = SVector(1.25)),
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0), SVector(dx)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(du)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    # One input notch between consecutive commands, starting and ending at rest.
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("bounded_input_variation"),
        OPDS.BoundedInputVariation(
            (u1, u2) -> maximum(abs.(u1 - u2)),
            du;
            target_input = SVector(0.0),
            initial_input = SVector(0.0),
        ),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    ctrl = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    @test ctrl !== nothing

    ctrl2 = roundtrip(ctrl)
    mem, mem2 = ST.initial_state(ctrl), ST.initial_state(ctrl2)
    @test mem2 == mem

    states = [SVector(x1) for x1 in (-1.875):dx:1.875]
    npts = 0
    for x in states
        @test ST.is_defined(ctrl2, mem2, x) == ST.is_defined(ctrl, mem, x)
        u = ST.output_control(ctrl, mem, x)
        u2 = ST.output_control(ctrl2, mem2, x)
        @test u == u2
        u === nothing || (npts += 1)
    end
    @test npts > 0

    # The guarantee that actually matters for deployment: a *fresh session*. An
    # in-process reload cannot fail on a closure-valued field, because the closure's
    # type is still defined in this session.
    ok, err, answers = fresh_session_answers(ctrl, states)
    ok || @info "fresh-session reload of the slew controller failed" err
    @test ok
    if ok
        @test length(answers) == length(states)
        for (i, x) in enumerate(states)
            @test answers[i] ==
                  (ST.is_defined(ctrl, mem, x), ST.output_control(ctrl, mem, x))
        end
    end
end

@testset "compact_controller keeps the answers and drops the transitions" begin
    # `SY.compact_controller` strips the abstraction's transition relation, which a
    # synthesized controller never reads again. It is what makes a controller small
    # enough to serialize, so it has to be answer-preserving.
    F_sys(x, u) = SVector(u[1], -0.5 * x[1])
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
    concrete_system = MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 2, 1, _X_, _U_)

    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        LazySets.Hyperrectangle(; low = SVector(-1.5, -0.5), high = SVector(-0.5, 0.5)),
        LazySets.Hyperrectangle(; low = SVector(0.5, -0.5), high = SVector(1.5, 0.5)),
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
    small = SY.compact_controller(ctrl)

    # The transition relation is gone ...
    @test SY.get_n_transitions(ST.domain(ctrl)) > 0
    @test SY.get_n_transitions(ST.domain(small)) == 0

    # ... and every answer is unchanged, before and after a round-trip.
    small2 = roundtrip(small)
    npts = 0
    for x1 in -1.9:0.4:1.9, x2 in -1.9:0.4:1.9
        x = SVector(x1, x2)
        @test ST.is_defined(small, nothing, x) == ST.is_defined(ctrl, nothing, x)
        @test ST.output_control(small, nothing, x) == ST.output_control(ctrl, nothing, x)
        @test ST.output_control(small2, nothing, x) == ST.output_control(ctrl, nothing, x)
        ST.output_control(ctrl, nothing, x) === nothing || (npts += 1)
    end
    @test npts > 0

    # A controller with no symbolic model to strip passes through untouched.
    affine = ST.AffineController(MS.AffineMap([1.0 0.0; 0.0 2.0], [0.5, -0.5]))
    @test SY.compact_controller(affine) === affine
end

end # module
