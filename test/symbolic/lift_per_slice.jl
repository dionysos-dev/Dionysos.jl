module TestLiftPerSlice

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS

# Base on a shared 2-state / 1-input grid with the given transitions (q′, q, u).
function build_base(transitions)
    Xmap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Xmap, (0,))   # q1 at x = 0.0
    MP.add_pos!(Xmap, (1,))   # q2 at x = 1.0
    Umap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Umap, (0,))
    base = SY.SymbolicModelList(Xmap, Umap)
    SY.add_transitions!(base.autom, transitions)
    return base
end

active_clock(tmax, tstep) = SY.ClockAbstraction(
    MS.ConstrainedLinearContinuousSystem([1.0;;], UT.box([0.0], [tmax])),
    tstep,
)

@testset "lift_per_slice - transitions follow each slice's dynamics" begin
    base_A = build_base([(2, 1, 1)])   # slice dynamics A: state 1 → 2
    base_B = build_base([(1, 2, 1)])   # slice dynamics B: state 2 → 1
    clock = active_clock(2.0, 1.0)     # tsteps [0, 1, 2] → 3 slices

    m = SY.lift_per_slice([base_A, base_B, base_A], clock)
    @test m isa SY.ClockLiftedSymbolicModel
    @test SY.get_n_state(m) == 2 * 3

    # Flattened ids for a few (x, t) coordinates.
    id(x, t) = SY.get_abstract_state(m, SVector(x, t))

    # Slice 1 uses base_A (1 → 2): expect (x=0,t=0) → (x=1,t=1).
    @test id(1.0, 1.0) in SY.post(m, id(0.0, 0.0), 1)
    # Slice 2 uses base_B (2 → 1): expect (x=1,t=1) → (x=0,t=2).
    @test id(0.0, 2.0) in SY.post(m, id(1.0, 1.0), 1)
    # Slice 2 does NOT have 1 → 2, so (x=0,t=1) has no successor (x=1,t=2).
    @test !(id(1.0, 2.0) in SY.post(m, id(0.0, 1.0), 1))
end

@testset "lift_per_slice - identical bases match lift(ClockLift)" begin
    base = build_base([(2, 1, 1), (2, 2, 1)])
    clock = active_clock(2.0, 1.0)

    m_slice = SY.lift_per_slice([base, base, base], clock)
    m_lift = SY.lift(SY.ClockLift(clock), base)

    @test SY.get_n_state(m_slice) == SY.get_n_state(m_lift)
    @test Set(SY.enum_transitions(m_slice)) == Set(SY.enum_transitions(m_lift))
end

end # module
