module TestMain

using Test
using StaticArrays
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

@testset "HierarchicalGridMapping" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    base = MP.ImplicitGridMapping(grid, (0, 0), (1, 1))  # 2x2 = 4 cells

    h = MP.HierarchicalGridMapping(base; nlevels = 2, div = 2)

    # It is a multi-level mapping, NOT a single-grid GridMapping.
    @test h isa MP.AbstractMultiLevelMapping
    @test h isa MP.AbstractMapping
    @test !(h isa MP.GridMapping)

    # Two levels: level 1 has 4 states, level 2 refines by 2 => 16 states.
    @test MP.get_n_state(h.levels[1]) == 4
    @test MP.get_n_state(h.levels[2]) == 16
    @test MP.get_n_state(h) == 20
    @test h.offsets == [0, 4]

    # Global labels of level 2 are offset past level 1.
    q2_local = 3
    q2_global = MP.get_state_by_pos(h, 2, MP.get_pos_by_state(h.levels[2], q2_local))
    @test q2_global == h.offsets[2] + q2_local

    # Global -> (level, pos) and Global -> coord agree with the located level.
    l, pos = MP.get_pos_by_state(h, q2_global)
    @test l == 2
    @test MP.get_coord_by_state(h, q2_global) ==
          MP.get_coord_by_state(h.levels[2], q2_local)

    # Set-based query on a level, offset into the global label space.
    rect = UT.box(SVector(-0.4, -0.4), SVector(0.4, 0.4))
    qs = MP.get_states_from_set(h, 1, rect, MP.OUTER)
    @test !isempty(qs)
    @test all(q -> q <= MP.get_n_state(h.levels[1]) + h.offsets[1], qs)
end

end # module
