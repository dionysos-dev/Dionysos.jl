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
grid_space = AB.NewGridSpaceHash(orig, h)
display(grid_space)

AB.add_to_gridspace_by_coords!(grid_space, (1.2, 3.5))
@test AB.get_gridspace_size(grid_space) == 1
AB.add_to_gridspace_by_coords!(grid_space, (-0.5002, -1.0))
@test AB.get_gridspace_size(grid_space) == 2
display(grid_space)

AB.add_to_gridspace!(grid_space, AB.HyperRectangle((1.0, 0.0), (11.0, 10.0)), AB.OUTER)
@test AB.get_gridspace_size(grid_space) == 67
println(grid_space.elems)

AB.remove_from_gridspace_by_coords!(grid_space, (2.0, 2.0))
@test AB.get_gridspace_size(grid_space) == 66
AB.remove_from_gridspace!(grid_space, AB.HyperRectangle((5.0, 5.0), (10000.0, 10000.0)), AB.INNER)
@test AB.get_gridspace_size(grid_space) == 48

pos_iter = AB.enumerate_gridspace_pos(grid_space)
@test length(pos_iter) == 48
display(pos_iter)

sub_set = AB.NewSubSet(grid_space)
AB.add_to_subset_all!(sub_set)
@test AB.get_subset_size(sub_set) == 48
AB.remove_from_subset!(sub_set, AB.HyperRectangle((1.0, 1.0), (2.0, 2.0)), AB.OUTER)
@test AB.get_subset_size(sub_set) == 45
AB.add_to_subset!(sub_set, AB.HyperRectangle((0.0, 0.0), (5.0, 5.0)), AB.INNER)
@test AB.get_subset_size(sub_set) == 46
display(sub_set)

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-2.0, 14.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.subset!(ax, 1:2, sub_set, fa = 0.3)
    Plot.set!(ax, 1:2, AB.HyperRectangle((1.0, 0.0), (8.0, 10.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((5.0, 5.0), (10000.0, 10000.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((1.0, 1.0), (2.0, 2.0)))
    Plot.set!(ax, 1:2, AB.HyperRectangle((0.0, 0.0), (5.0, 5.0)))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
