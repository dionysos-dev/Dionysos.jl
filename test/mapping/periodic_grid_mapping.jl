module TestPeriodicGridMapping

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# Helper: build a grid with correct periodic origin alignment
function _aligned_grid(periodic_dims, periods, start, h, N::Int)
    orig = zeros(Float64, N)
    for (i, d) in enumerate(periodic_dims)
        orig[d] = start[i] + h[d] / 2
    end
    # non-periodic dims stay at 0.0
    return MP.GridFree(SVector{N}(orig), h)
end

@testset "PeriodicGridMapping metadata + wrapping (2D, dim 1 periodic)" begin
    periodic_dims = SVector(1)        # periodic in x
    periods = SVector(4.0)
    start = SVector(0.0)
    h = SVector(1.0, 1.0)

    grid = _aligned_grid(periodic_dims, periods, start, h, 2)

    # underlying explicit mapping so we can test add_pos!/add_set!
    base = MP.ExplicitGridMapping(grid)
    pm = MP.PeriodicGridMapping(periodic_dims, periods, start, base)

    @test MP.get_periodic_dims(pm) == periodic_dims
    @test MP.get_periods(pm) == periods
    @test MP.get_periodic_starts(pm) == start

    @test MP.is_periodic(pm)
    @test MP.is_periodic(pm, 1)
    @test !MP.is_periodic(pm, 2)

    # periodic index map sanity
    @test pm.periodic_index_map == (1, nothing)

    # wrap_pos (span = period/h = 4/1 = 4)
    @test MP.wrap_pos(pm, (0, 7)) == (0, 7)
    @test MP.wrap_pos(pm, (4, 7)) == (0, 7)
    @test MP.wrap_pos(pm, (5, 7)) == (1, 7)
    @test MP.wrap_pos(pm, (-1, 7)) == (3, 7)

    # wrap_coord
    x = SVector(4.1, 2.2)
    wx = MP.wrap_coord(pm, x)
    @test isapprox(wx[1], 0.1; atol = 1e-9)
    @test isapprox(wx[2], 2.2; atol = 1e-12)

    xneg = SVector(-1.0, 0.0)
    wxneg = MP.wrap_coord(pm, xneg)
    @test isapprox(wxneg[1], 3.0; atol = 1e-9)
end

@testset "PeriodicGridMapping delegates get_state_by_pos / is_valid_pos / add_pos!" begin
    periodic_dims = SVector(1)
    periods = SVector(4.0)
    start = SVector(0.0)
    h = SVector(1.0, 1.0)
    grid = _aligned_grid(periodic_dims, periods, start, h, 2)

    base = MP.ExplicitGridMapping(grid)

    # Fill base with a few positions so states exist
    q00 = MP.add_pos!(base, (0, 0))
    q10 = MP.add_pos!(base, (1, 0))
    q30 = MP.add_pos!(base, (3, 0))

    pm = MP.PeriodicGridMapping(periodic_dims, periods, start, base)

    # validity wraps before checking underlying
    @test MP.is_valid_pos(pm, (0, 0))
    @test MP.is_valid_pos(pm, (4, 0))     # wraps to (0,0) which is valid
    @test MP.is_valid_pos(pm, (-1, 0))    # wraps to (3,0) which is valid

    # get_state_by_pos wraps
    @test MP.get_state_by_pos(pm, (0, 0)) == q00
    @test MP.get_state_by_pos(pm, (4, 0)) == q00
    @test MP.get_state_by_pos(pm, (5, 0)) == q10
    @test MP.get_state_by_pos(pm, (-1, 0)) == q30

    # add_pos! wraps then inserts into underlying
    q = MP.add_pos!(pm, (4, 2))  # wraps to (0,2)
    @test MP.get_state_by_pos(pm, (0, 2)) == q
    @test MP.get_state_by_pos(pm, (4, 2)) == q
end

@testset "PeriodicGridMapping coord<->state path" begin
    periodic_dims = SVector(1)
    periods = SVector(4.0)
    start = SVector(0.0)
    h = SVector(1.0, 1.0)
    grid = _aligned_grid(periodic_dims, periods, start, h, 2)

    base = MP.ExplicitGridMapping(grid)
    q = MP.add_pos!(base, (0, 0))
    pm = MP.PeriodicGridMapping(periodic_dims, periods, start, base)

    # coordinate that is outside period in dim 1 should wrap
    # With your GridFree convention, coord-by-pos for pos=(0,0) should be orig + pos*h
    # orig[1] = start + h/2 = 0.5, so coord is (0.5, 0.0)
    x = SVector(4.5, 0.0)  # wraps to 0.5 in dim1
    @test MP.get_state_by_coord(pm, x) == q

    # coord-by-state is computed via pos-by-state + grid coord
    xq = MP.get_coord_by_state(pm, q)
    @test isapprox(xq[1], 0.5; atol = 1e-9)
end

@testset "PeriodicGridMapping add_set! (Explicit underlying only)" begin
    periodic_dims = SVector(1)
    periods = SVector(4.0)
    start = SVector(0.0)
    h = SVector(1.0, 1.0)
    grid = _aligned_grid(periodic_dims, periods, start, h, 2)

    base = MP.ExplicitGridMapping(grid)
    pm = MP.PeriodicGridMapping(periodic_dims, periods, start, base)

    rect = UT.box(SVector(0.0, 0.0), SVector(2.0, 0.0))
    MP.cover!(pm, rect, MP.OUTER)

    # should have inserted positions (0,0), (1,0), (2,0) (mod 4 in dim 1)
    @test MP.is_valid_pos(pm, (0, 0))
    @test MP.is_valid_pos(pm, (1, 0))
    @test MP.is_valid_pos(pm, (2, 0))

    # wrapped equivalents should also be valid if they wrap to inserted positions
    @test MP.is_valid_pos(pm, (4, 0))  # wraps to (0,0)
end

@testset "PeriodicGridMapping constructor checks (origin alignment / divisibility)" begin
    periodic_dims = SVector(1)
    periods = SVector(4.0)
    start = SVector(0.0)

    # bad origin (should be start + h/2 = 0.5 in dim1)
    bad_grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    bad_base = MP.ExplicitGridMapping(bad_grid)
    @test_throws ErrorException MP.PeriodicGridMapping(
        periodic_dims,
        periods,
        start,
        bad_base,
    )

    # bad divisibility: period/h not integer
    hbad = SVector(3.0, 1.0)
    grid_bad_step = _aligned_grid(periodic_dims, periods, start, hbad, 2)
    base_bad_step = MP.ExplicitGridMapping(grid_bad_step)
    @test_throws ErrorException MP.PeriodicGridMapping(
        periodic_dims,
        periods,
        start,
        base_bad_step,
    )
end

end # module
