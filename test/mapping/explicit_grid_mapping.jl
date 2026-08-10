module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots
import LazySets

@testset "ExplicitGridMapping" begin
    orig = SVector(0.0, 0.0)
    h = SVector(1.0, 2.0)

    grid = MP.GridFree(orig, h)

    # --- basic grid conversion sanity ---
    @test MP.get_origin(grid) == orig
    @test MP.get_h(grid) == h

    x = SVector(1.2, 3.5)
    pos = MP.get_pos_by_coord(grid, x)
    @test pos == (1, 2)  # round(1.2/1)=1, round(3.5/2)=2

    x2 = MP.get_coord_by_pos(grid, pos)
    @test x2 == SVector(1.0, 4.0)

    # --- explicit mapping starts empty ---
    m = MP.ExplicitGridMapping(grid)
    @test MP.get_n_state(m) == 0

    # add two positions (dedup check)
    q1 = MP.add_pos!(m, (1, 2))
    @test q1 == 1
    @test MP.get_n_state(m) == 1
    @test MP.get_pos_by_state(m, 1) == (1, 2)
    @test MP.get_state_by_pos(m, (1, 2)) == 1

    q1b = MP.add_pos!(m, (1, 2))
    @test q1b == 1
    @test MP.get_n_state(m) == 1

    q2 = MP.add_pos!(m, (-1, -1))
    @test q2 == 2
    @test MP.get_n_state(m) == 2

    # --- cover! using OUTER inclusion ---
    rect = LazySets.Hyperrectangle(; low = SVector(1.0, 0.0), high = SVector(11.0, 10.0))
    MP.cover!(m, rect, MP.OUTER)

    @test MP.get_n_state(m) == 67

    # --- plotting smoke test (if you kept Mapping recipes) ---
    fig = plot(; aspect_ratio = :equal)
    plot!(fig, m)
    @test isa(fig, Plots.Plot)

    m2 = MP.ExplicitGridMapping(MP.GridFree(orig, SVector(0.1, 0.1)))
    # compile once
    @allocated MP.cover!(m2, rect, MP.OUTER)

    alloc2 = @allocated MP.cover!(m2, rect, MP.OUTER)
    @test alloc2 == 0
end

@testset "get_states_from_set_strict (allin flag over a partial domain)" begin
    # 1×1 cells centered at integer coords; map exactly the 2×2 block covering
    # [0, 1]²  → cells (0,0), (1,0), (0,1), (1,1).
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    MP.cover!(
        m,
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(1.0, 1.0)),
        MP.OUTER,
    )
    n = MP.get_n_state(m)
    @test n == 4

    # a set inside the mapped domain: every covered cell is a valid state
    inside = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(1.0, 1.0))
    qs, allin = MP.get_states_from_set_strict(m, inside, MP.OUTER)
    @test allin
    @test Set(qs) == Set(1:n)

    # a set spilling outside the mapped domain: allin drops to false, and only
    # the mapped cells become states
    outside = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(5.0, 5.0))
    qs2, allin2 = MP.get_states_from_set_strict(m, outside, MP.OUTER)
    @test !allin2
    @test Set(qs2) == Set(1:n)
    # the non-strict form is just the states, flag dropped
    @test Set(MP.get_states_from_set(m, outside, MP.OUTER)) == Set(qs2)

    # EmptySet → no states, vacuously all-in
    qe, aine = MP.get_states_from_set_strict(m, LazySets.EmptySet(2), MP.OUTER)
    @test isempty(qe) && aine

    # single coordinate: inside → ([q], true); outside → (nothing, false)
    qc, ainc = MP.get_states_from_set_strict(m, SVector(0.0, 0.0), MP.OUTER)
    @test ainc && qc == [MP.get_state_by_pos(m, (0, 0))]
    qc2, ainc2 = MP.get_states_from_set_strict(m, SVector(9.0, 9.0), MP.OUTER)
    @test qc2 === nothing && !ainc2

    # a union routes through the per-subset strict form and deduplicates
    u = UT.set_union([
        LazySets.Hyperrectangle(; low = SVector(-0.1, -0.1), high = SVector(0.1, 0.1)),  # cell (0,0)
        LazySets.Hyperrectangle(; low = SVector(0.9, 0.9), high = SVector(1.1, 1.1)),    # cell (1,1)
    ])
    qu, ainu = MP.get_states_from_set_strict(m, u, MP.OUTER)
    @test ainu
    @test Set(qu) == Set([MP.get_state_by_pos(m, (0, 0)), MP.get_state_by_pos(m, (1, 1))])
    @test length(qu) == length(unique(qu))

    # a set_minus carves the hole with the inverted inclusion mode
    sm = UT.set_minus(
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(1.0, 1.0)),      # kept: all 4 cells (OUTER)
        LazySets.Hyperrectangle(; low = SVector(0.4, 0.4), high = SVector(1.6, 1.6)),      # hole: cell (1,1) (INNER)
    )
    qm, _ = MP.get_states_from_set_strict(m, sm, MP.OUTER)
    @test MP.get_state_by_pos(m, (1, 1)) ∉ qm    # carved out
    @test MP.get_state_by_pos(m, (0, 0)) ∈ qm    # kept
end

end # module TestMain
