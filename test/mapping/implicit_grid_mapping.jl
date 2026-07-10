module TestImplicitGridMapping

using Test
using StaticArrays
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

@testset "ImplicitGridMapping basic indexing" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))

    min_pos = (0, 0)
    max_pos = (2, 3)  # 3 x 4 box
    m = MP.ImplicitGridMapping(grid, min_pos, max_pos)

    @test MP.get_grid(m) === grid
    @test m.n_per_dim == (3, 4)
    @test m.strides == (1, 3)   # row-major: x stride 1, y stride n_x = 3
    @test MP.get_n_state(m) == 12
    @test collect(MP.enum_states(m)) == collect(1:12)

    # --- get_state_by_pos matches row-major expectation ---
    # q = 1 + (x-minx)*stride1 + (y-miny)*stride2
    @test MP.get_state_by_pos(m, (0, 0)) == 1
    @test MP.get_state_by_pos(m, (1, 0)) == 2
    @test MP.get_state_by_pos(m, (2, 0)) == 3
    @test MP.get_state_by_pos(m, (0, 1)) == 4
    @test MP.get_state_by_pos(m, (2, 1)) == 6
    @test MP.get_state_by_pos(m, (0, 3)) == 10
    @test MP.get_state_by_pos(m, (2, 3)) == 12

    # --- inverse mapping ---
    @test MP.get_pos_by_state(m, 1) == (0, 0)
    @test MP.get_pos_by_state(m, 2) == (1, 0)
    @test MP.get_pos_by_state(m, 3) == (2, 0)
    @test MP.get_pos_by_state(m, 4) == (0, 1)
    @test MP.get_pos_by_state(m, 12) == (2, 3)

    # --- round-trip for all states ---
    for q in MP.enum_states(m)
        pos = MP.get_pos_by_state(m, q)
        @test MP.get_state_by_pos(m, pos) == q
    end

    # --- round-trip for all positions ---
    for x in 0:2, y in 0:3
        pos = (x, y)
        q = MP.get_state_by_pos(m, pos)
        @test MP.get_pos_by_state(m, q) == pos
    end

    # --- out-of-range checks ---
    @test MP.get_state_by_pos(m, (-1, 0)) == nothing
    @test MP.get_state_by_pos(m, (3, 0)) == nothing
    @test MP.get_state_by_pos(m, (0, 4)) == nothing
    @test MP.get_pos_by_state(m, 0) == nothing
end

@testset "ImplicitGridMapping from rectangle" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 2.0))

    # rectangle in REAL coords
    rect = UT.box(SVector(1.0, 0.0), SVector(11.0, 10.0))

    m = MP.ImplicitGridMapping(grid, rect; incl_mode = MP.OUTER)

    # Expected OUTER pos limits for your get_pos_lims_outer:
    # x in [1,11] with h=1 => pos_x 1..11  (11)
    # y in [0,10] with h=2 => pos_y 0..5   (6)
    @test m.min_pos == (1, 0)
    @test m.max_pos == (11, 5)
    @test m.n_per_dim == (11, 6)
    @test MP.get_n_state(m) == 66

    # sanity: corners
    @test MP.get_state_by_pos(m, (1, 0)) == 1
    @test MP.get_state_by_pos(m, (11, 5)) == 66
    @test MP.get_pos_by_state(m, 66) == (11, 5)
end

end # module
