include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "GridSpaceSubSet" begin
orig = (0.0, 0.0)
h = (1.0, 2.0)
gridspace = AB.NewGridSpaceList(orig, h)

AB.add_coord!(gridspace, (1.2, 3.5))
@test AB.get_ncells(gridspace) == 1
AB.add_coord!(gridspace, (-0.5002, -1.0))
@test AB.get_ncells(gridspace) == 2

AB.add_set!(gridspace, AB.HyperRectangle((1.0, 0.0), (11.0, 10.0)), AB.OUTER)
@test AB.get_ncells(gridspace) == 67

AB.remove_coord!(gridspace, (2.0, 2.0))
@test AB.get_ncells(gridspace) == 66
AB.remove_set!(gridspace, AB.HyperRectangle((5.0, 5.0), (10000.0, 10000.0)), AB.INNER)
@test AB.get_ncells(gridspace) == 48

pos_iter = AB.enum_pos(gridspace)
@test length(pos_iter) == 48

subset = AB.NewSubSet(gridspace)
AB.add_all!(subset)
@test AB.get_ncells(subset) == 48
AB.remove_set!(subset, AB.HyperRectangle((1.0, 1.0), (2.0, 2.0)), AB.OUTER)
@test AB.get_ncells(subset) == 45
AB.add_set!(subset, AB.HyperRectangle((0.0, 0.0), (5.0, 5.0)), AB.INNER)
@test AB.get_ncells(subset) == 46

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-2.0, 14.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.subset!(ax, 1:2, subset, fa = 0.3)
    Plot.set!(ax, 1:2, AB.HyperRectangle((1.0, 0.0), (8.0, 10.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((5.0, 5.0), (10000.0, 10000.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((1.0, 1.0), (2.0, 2.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((0.0, 0.0), (5.0, 5.0)))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
