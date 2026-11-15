module TestMain
using Test
using StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

sleep(0.1) # used for good printing
println("Started test")

@testset "PeriodicDomainList Construction (dims 1 and 3 periodic)" begin
    # 3D setting
    periodic_dims = SVector(1, 3)
    start = SVector(0.0, 0.0)
    periods = SVector(4.0, 8.0)
    h = SVector(1.0, 1.0, 2.0)

    # Grid origin: ensure orig[d] = start[i] + h[d]/2 for periodic dims
    orig = Vector{Float64}(undef, 3)
    orig[1] = start[1] + h[1] / 2  # for dim 1
    orig[2] = 0.0                  # non-periodic dim 2
    orig[3] = start[2] + h[3] / 2  # for dim 3
    grid1 = DO.GridFree(SVector(orig...), h)

    domain1 = DO.PeriodicDomainList(periodic_dims, periods, start, grid1)

    @test DO.get_periodic_dims(domain1) == periodic_dims
    @test DO.get_periods(domain1) == periods
    @test DO.get_periodic_starts(domain1) == start
    @test DO.is_periodic(domain1)
    @test DO.is_periodic(domain1, 1)
    @test !DO.is_periodic(domain1, 2)
    @test DO.is_periodic(domain1, 3)

    # Wrapping behavior
    @test DO.wrap_pos(domain1, (4, 8, 8)) == (0, 8, 0)
    @test DO.wrap_pos(domain1, (2, 5, 10)) == (2, 5, 2)

    # --- Constructor: from (start, h) ---
    domain2 = DO.PeriodicDomainList(periodic_dims, periods, start, h)
    @test DO.has_same_periodicity(domain1, domain2)
    @test DO.wrap_pos(domain2, (4, 2, 8)) == (0, 2, 0)
    @test DO.wrap_pos(domain2, (5, 2, 10)) == (1, 2, 2)

    # --- Constructor: from (periods, h) — assumes start = 0 ---
    domain3 = DO.PeriodicDomainList(periodic_dims, periods, h)
    @test DO.has_same_periodicity(domain1, domain3)
    @test DO.get_periodic_starts(domain3) == SVector(0.0, 0.0)
    @test DO.wrap_pos(domain3, (8, 1, 13)) == (0, 1, 1)

    # --- Constructor: non-periodic from h ---
    domain4 = DO.PeriodicDomainList(h)
    @test !DO.has_same_periodicity(domain1, domain4)
    @test DO.get_periodic_dims(domain4) == SVector{0, Int}()
    @test !DO.is_periodic(domain4)

    # --- Constructor: non-periodic from grid ---
    grid2 = DO.GridFree(SVector(0.0, 0.0, 0.0), h)
    domain5 = DO.PeriodicDomainList(grid2)
    @test !DO.has_same_periodicity(domain1, domain5)
    @test DO.has_same_periodicity(domain4, domain5)
    @test !DO.is_periodic(domain5)
    @test DO.get_periodic_dims(domain5) == SVector{0, Int}()
end

@testset "wrap_pos and wrap_coord with periodic dims 1, 3, 5 (including negatives)" begin
    # 5D setup
    periodic_dims = SVector(1, 3, 5)
    start = SVector(0.0, 0.0, 0.0)
    periods = SVector(4.0, 6.0, 10.0)
    h = SVector(1.0, 1.0, 1.0, 1.0, 2.0)  # note: dim 5 has larger step

    # Construct grid origin to satisfy periodic alignment
    orig = Vector{Float64}(undef, 5)
    for (i, d) in enumerate(periodic_dims)
        orig[d] = start[i] + h[d] / 2
    end
    for d in setdiff(1:5, periodic_dims)
        orig[d] = 0.0
    end
    grid = DO.GridFree(SVector(orig...), h)
    domain = DO.PeriodicDomainList(periodic_dims, periods, start, grid)

    # wrap_pos tests (NTuple{5, Int})
    @test DO.wrap_pos(domain, (4, 2, 6, 1, 10)) == (0, 2, 0, 1, 0)
    @test DO.wrap_pos(domain, (-1, 0, -1, 0, -1)) == (3, 0, 5, 0, 4)
    @test DO.wrap_pos(domain, (9, 1, 14, 3, 22)) == (1, 1, 2, 3, 2)

    # wrap_coord tests (SVector{5, Float64})
    coord = SVector(4.1, 2.2, 6.3, 0.5, 10.5)
    wrapped = DO.wrap_coord(domain, coord)
    @test isapprox(wrapped[1], 0.1; atol = 1e-8)
    @test wrapped[2] ≈ 2.2
    @test isapprox(wrapped[3], 0.3; atol = 1e-8)
    @test wrapped[4] ≈ 0.5
    @test isapprox(wrapped[5], 0.5; atol = 1e-8)

    coord_neg = SVector(-1.0, 1.0, -2.5, 1.0, -5.0)
    wrapped_neg = DO.wrap_coord(domain, coord_neg)
    @test isapprox(wrapped_neg[1], 3.0; atol = 1e-8)
    @test wrapped_neg[2] ≈ 1.0
    @test isapprox(wrapped_neg[3], 3.5; atol = 1e-8)
    @test wrapped_neg[4] ≈ 1.0
    @test isapprox(wrapped_neg[5], 5.0; atol = 1e-8)
end

@testset "Membership and set logic (periodic dims 1 and 3)" begin
    # Domain settings
    periodic_dims = SVector(1, 3)
    start = SVector(0.0, 0.0)
    periods = SVector(4.0, 8.0)
    h = SVector(1.0, 1.0, 2.0)

    # Grid origin aligned with periodicity
    orig = Vector{Float64}(undef, 3)
    orig[1] = start[1] + h[1] / 2
    orig[2] = 0.0
    orig[3] = start[2] + h[3] / 2
    grid = DO.GridFree(SVector(orig...), h)

    domain = DO.PeriodicDomainList(periodic_dims, periods, start, grid)

    # Add a wrapped position (gets wrapped from (4, 1, 8) → (0, 1, 0))
    DO.add_pos!(domain, (4, 1, 8))
    @test (0, 1, 0) in domain
    @test (4, 1, 8) in domain  # Should still pass due to wrapping
    @test (8, 1, 16) in domain  # Multiple wraps
    @test !((2, 3, 4) in domain)  # Never added

    # Add a negative-wrapped position
    DO.add_pos!(domain, (-1, 1, -2))  # wraps to (3, 1, 6)
    @test (3, 1, 6) in domain
    @test (-1, 1, -2) in domain

    # Add a wrapped position
    DO.add_pos!(domain, (4, 1, 8))  # wraps to (0, 1, 0)

    @test (0, 1, 0) in domain
    @test (4, 1, 8) in domain

    # Remove using wrapped form
    DO.remove_pos!(domain, (0, 1, 0))
    @test !((0, 1, 0) in domain)
    @test !((4, 1, 8) in domain)

    # Add again
    DO.add_pos!(domain, (8, 1, 16))  # wraps to same (0, 1, 0)
    @test (8, 1, 16) in domain

    # Remove using unwrapped form
    DO.remove_pos!(domain, (8, 1, 16))
    @test !((0, 1, 0) in domain)
    @test !((8, 1, 16) in domain)
end

@testset "Set operations (periodic dims 1 and 3)" begin
    # Settings
    periodic_dims = SVector(1, 3)
    start = SVector(0.0, 0.0)
    periods = SVector(4.0, 8.0)
    h = SVector(1.0, 1.0, 2.0)

    # Grid aligned with periodic dims
    orig = Vector{Float64}(undef, 3)
    orig[1] = start[1] + h[1]/2
    orig[2] = 0.0
    orig[3] = start[2] + h[3]/2
    grid = DO.GridFree(SVector(orig...), h)

    d1 = DO.PeriodicDomainList(periodic_dims, periods, start, grid)
    d2 = DO.PeriodicDomainList(periodic_dims, periods, start, grid)

    # (4, 1, 8) wraps to (0, 1, 0)
    DO.add_pos!(d1, (0, 1, 0))
    DO.add_pos!(d2, (4, 1, 8))

    @test DO.issubset(d2, d1)

    DO.union!(d1, d2)
    @test (0, 1, 0) in d1
    @test (4, 1, 8) in d1  # because it's the same point wrapped

    DO.setdiff!(d1, d2)

    @test !((0, 1, 0) in d1)
    @test !((4, 1, 8) in d1)
end

@testset "Interface Methods Delegated or Wrapped (periodic dims 1 and 3)" begin

    # Setup: 3D grid with periodic dims 1 and 3
    periodic_dims = SVector(1, 3)
    start = SVector(0.0, 0.0)
    periods = SVector(4.0, 8.0)
    h = SVector(1.0, 1.0, 2.0)

    orig = Vector{Float64}(undef, 3)
    orig[1] = start[1] + h[1] / 2
    orig[2] = 0.0
    orig[3] = start[2] + h[3] / 2
    grid = DO.GridFree(SVector(orig...), h)

    domain = DO.PeriodicDomainList(periodic_dims, periods, start, grid)

    # Fill with some points
    DO.add_pos!(domain, (4, 0, 8))  # wraps to (0, 0, 0)
    DO.add_pos!(domain, (1, 1, 2))

    # get_grid
    @test DO.get_grid(domain) === grid

    # enum_pos
    positions = collect(DO.enum_pos(domain))
    @test (0, 0, 0) in positions
    @test (1, 1, 2) in positions
    @test length(positions) == 2

    # get_ncells
    @test DO.get_ncells(domain) == 2

    # get_somepos
    somepos = DO.get_somepos(domain)
    @test somepos in positions

    # isempty
    @test !isempty(domain)
    DO.empty!(domain)
    @test isempty(domain)

    # in(pos, domain) with wrapping
    DO.add_pos!(domain, (4, 0, 8))  # wrapped
    @test (0, 0, 0) in domain
    @test (4, 0, 8) in domain

    # get_pos_by_coord (coord → pos → wrapped)
    coord = SVector(4.5, 0.0, 8.5)
    pos = DO.get_pos_by_coord(domain, coord)
    @test pos == (0, 0, 0)

    # get_coord_by_pos (pos → wrapped → coord)
    coord_back = DO.get_coord_by_pos(domain, (4, 0, 8))
    @test isapprox(coord_back[1], 0.5; atol = 1e-8)
    @test isapprox(coord_back[3], 1.0; atol = 1e-8)

    # get_subset_pos: bounding box around wrapped content
    INCL_MODE = DO.OUTER
    rect = UT.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    subset = DO.get_subset_pos(domain, rect, INCL_MODE)
    @test (0, 0, 0) in subset
end

@testset "Periodic set operations (add/remove/subset)" begin
    # Domain setup
    periodic_dims = SVector(1, 3)
    start = SVector(0.0, 0.0)
    periods = SVector(4.0, 8.0)
    h = SVector(1.0, 1.0, 2.0)

    # Aligned grid
    origin = SVector(start[1] + h[1]/2, 0.0, start[2] + h[3]/2)
    grid = DO.GridFree(origin, h)
    domain = DO.PeriodicDomainList(periodic_dims, periods, start, grid)

    # ---- Test: add_set! with a HyperRectangle
    rect = UT.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(2.0, 1.0, 2.0))
    DO.add_set!(domain, rect, DO.OUTER)
    @test (0, 0, 0) in domain
    @test (1, 0, 0) in domain
    @test (2, 0, 0) in domain

    # ---- Test: remove_set! removes subset
    DO.remove_set!(domain, rect, DO.OUTER)
    @test !((0, 0, 0) in domain)
    @test !((1, 0, 0) in domain)

    # ---- Test: remove_coord!
    pt = SVector(1.5, 0.0, 0.0)
    DO.add_set!(domain, rect, DO.OUTER)
    @test (1, 0, 0) in domain
    DO.remove_coord!(domain, pt)
    @test !((1, 0, 0) in domain)

    # ---- Test: add_set! with LazyUnionSetArray
    r1 = UT.HyperRectangle(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    r2 = UT.HyperRectangle(SVector(1.0, -1.0, 1.0), SVector(3.0, 1.0, 3.0))
    union = UT.LazyUnionSetArray([r1, r2])
    DO.add_set!(domain, union, DO.OUTER)
    @test (0, 0, 0) in domain
end

sleep(0.1) # used for good printing
println("End test")
end # end module
