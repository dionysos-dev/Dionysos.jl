module TestMain

using Test
using StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

sleep(0.1)
println("Started test")

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
    rect = UT.box(SVector(1.0, 0.0), SVector(11.0, 10.0))
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

sleep(0.1)
println("End test")

end # module TestMain
