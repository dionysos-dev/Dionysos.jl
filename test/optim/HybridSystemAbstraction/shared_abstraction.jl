module TestSharedAbstraction

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI
import MathematicalSystems as MS
using HybridSystems
import LazySets

const HSA = AB.HybridSystemAbstraction

const XBOX = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
const UBOX = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))

_mode(f, X = XBOX) = MS.ConstrainedBlackBoxControlContinuousSystem(f, 1, 1, X, UBOX)

# Two dynamics that compute the same thing but are *not* `===`: a system is an
# immutable struct, so `===` compares structurally, and two closures of different
# types never match. This is the case a declaration is for — the user knows the
# two modes are equivalent, the framework cannot see it.
const F_PLAIN = let a = 1.0
    (x, u) -> a .* u
end
const F_EQUIVALENT = let a = 1.0, b = 0.0
    (x, u) -> a .* u .+ b
end

function _two_mode_system(m1, m2)
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    HybridSystems.add_transition!(automaton, 2, 1, 2)
    guard = LazySets.Hyperrectangle(; low = SVector(0.5), high = SVector(1.0))
    return HybridSystems.HybridSystem(
        automaton,
        [m1, m2],
        [ST.GuardedResetMap(guard), ST.GuardedResetMap(guard)],
        [HybridSystems.AutonomousSwitching(), HybridSystems.AutonomousSwitching()],
    )
end

# Counts how many times a fresh uniform-grid abstraction is actually built
# (`optimizer_list` is typed `Vector{Function}`, hence a closure, not a struct).
function counting_factory()
    calls = Ref(0)
    factory = function ()
        calls[] += 1
        return MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    end
    return factory, calls
end

const KWARGS = Dict(
    "state_grid" => MP.GridFree(SVector(0.0), SVector(0.2)),
    "input_grid" => MP.GridFree(SVector(0.0), SVector(0.5)),
    "time_step" => 0.1,
    "approx_mode" => AB.UniformGridAbstraction.CENTER_SIMULATION,
    "print_level" => 0,
)

function _build(hs; shared_abstraction = nothing, kwargs = [copy(KWARGS), copy(KWARGS)])
    factory, calls = counting_factory()
    models = HSA.build_mode_symbolic_models(
        hs,
        Function[factory, factory],
        kwargs;
        shared_abstraction = shared_abstraction,
    )
    return models, calls[]
end

@testset "declared sharing builds the abstraction once" begin
    hs = _two_mode_system(_mode(F_PLAIN), _mode(F_EQUIVALENT))

    models_plain, calls_plain = _build(hs)
    @test calls_plain == 2

    models_shared, calls_shared = _build(hs; shared_abstraction = [nothing, 1])
    @test calls_shared == 1
    # The shared mode is the very same model object, so nothing was recomputed.
    @test models_shared[2] === models_shared[1]

    # ... and the abstraction is the one a plain build would have produced.
    @test SY.get_n_state(models_shared[2]) == SY.get_n_state(models_plain[2])
    @test length(collect(SY.enum_transitions(models_shared[2]))) ==
          length(collect(SY.enum_transitions(models_plain[2])))
end

@testset "identical modes are reused without a declaration" begin
    # A system is an immutable struct, so `===` is structural: whether the two
    # modes are the same object or two identically-built ones, reuse is provably
    # safe and happens automatically. This is the walking model's case.
    shared_object = _mode(F_PLAIN)
    _, calls_same_object = _build(_two_mode_system(shared_object, shared_object))
    @test calls_same_object == 1

    _, calls_rebuilt = _build(_two_mode_system(_mode(F_PLAIN), _mode(F_PLAIN)))
    @test calls_rebuilt == 1
end

@testset "invalid sharing declarations are refused" begin
    hs = _two_mode_system(_mode(F_PLAIN), _mode(F_EQUIVALENT))

    # Forward reference: a mode may only reuse an earlier one.
    @test_throws ErrorException _build(hs; shared_abstraction = [2, nothing])
    # Self reference.
    @test_throws ErrorException _build(hs; shared_abstraction = [nothing, 2])
    # Wrong length.
    @test_throws ErrorException _build(hs; shared_abstraction = [nothing])

    # Different optimizer configuration on the two modes.
    other = copy(KWARGS)
    other["state_grid"] = MP.GridFree(SVector(0.0), SVector(0.1))
    @test_throws ErrorException _build(
        hs;
        shared_abstraction = [nothing, 1],
        kwargs = [copy(KWARGS), other],
    )
end

@testset "modes with different state sets cannot share" begin
    wider = LazySets.Hyperrectangle(; low = SVector(-2.0), high = SVector(2.0))
    hs = _two_mode_system(_mode(F_PLAIN), _mode(F_PLAIN, wider))

    # Genuinely different domains: abstracted separately, and a declaration that
    # they are shareable is refused rather than silently honoured.
    _, calls = _build(hs)
    @test calls == 2
    @test_throws ErrorException _build(hs; shared_abstraction = [nothing, 1])
end

end # module
