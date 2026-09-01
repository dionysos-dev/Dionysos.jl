module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import HybridSystems
import LazySets

@testset "environment_input reads the ecosystem's ownership declarations" begin
    A1 = [0.5 0.0; 0.0 0.4]
    A2 = [0.3 0.1; 0.0 0.2]
    X = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [1.0, 1.0])
    U = LazySets.Hyperrectangle(; low = [-0.5], high = [0.5])
    W = LazySets.Hyperrectangle(; low = [-0.1, -0.1], high = [0.1, 0.1])

    # No noise declared: the controller owns everything, whatever else the type carries.
    plain = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        (x, u) -> A1 * x,
        2,
        1,
        X,
        U,
    )
    @test ST.environment_input(plain) === nothing

    # A declared noise set is the environment's choice set, verbatim.
    noisy = MathematicalSystems.NoisyConstrainedLinearControlDiscreteSystem(
        A1,
        reshape([0.0, 1.0], 2, 1),
        [1.0 0.0; 0.0 1.0],
        X,
        U,
        W,
    )
    @test ST.environment_input(noisy) === W

    # `discreteswitchedsystem` declares autonomous switching, under which the modes are the
    # environment's; making the switching controlled hands them back.
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    @test ST.environment_input(f) == 1:2

    f_controlled = ST.with_switching(f, HybridSystems.ControlledSwitching())
    @test ST.environment_input(f_controlled) === nothing
    # `with_switching` re-declares; it does not rebuild.
    @test f_controlled.automaton === f.automaton
    @test f_controlled.resetmaps === f.resetmaps

    # Mixed ownership has no whole-system answer and must be refused, not guessed.
    mixed = HybridSystems.HybridSystem(
        f.automaton,
        f.modes,
        f.resetmaps,
        [HybridSystems.AutonomousSwitching(), HybridSystems.ControlledSwitching()],
        f.ext,
    )
    @test_throws ErrorException ST.environment_input(mixed)

    # Not a system type at all: the default is pure synthesis, silently.
    @test ST.environment_input(:not_a_system) === nothing
end

end # module TestMain
