include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays

using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "GridSpace" begin
orig = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid_space = AB.NewGridSpaceHash(orig, h)
display(grid_space)

AB.add_to_gridspace_by_coords!(grid_space, SVector(1.2, 3.5))
@test AB.get_gridspace_size(grid_space) == 1
AB.add_to_gridspace_by_coords!(grid_space, SVector(-0.5002, -1.0))
@test AB.get_gridspace_size(grid_space) == 2
display(grid_space)

AB.add_to_gridspace!(grid_space, AB.HyperRectangle(SVector(1.0, 0.0), SVector(11.0, 10.0)), AB.OUTER)
@test AB.get_gridspace_size(grid_space) == 67
println(grid_space.elems)

AB.remove_from_gridspace_by_coords!(grid_space, SVector(2.0, 2.0))
@test AB.get_gridspace_size(grid_space) == 66
AB.remove_from_gridspace!(grid_space, AB.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)), AB.INNER)
@test AB.get_gridspace_size(grid_space) == 48
println(grid_space.elems)

pos_coll = AB.enumerate_gridspace_pos(grid_space)
@test length(pos_coll) == 48
display(pos_coll)

sub_space = AB.NewSubSpace(grid_space)
AB.add_to_subspace_all!(sub_space)
@test AB.get_subspace_size(sub_space) == 48
AB.remove_from_subspace!(sub_space, AB.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)), AB.OUTER)
@test AB.get_subspace_size(sub_space) == 45
AB.add_to_subspace!(sub_space, AB.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)), AB.INNER)
@test AB.get_subspace_size(sub_space) == 46
display(sub_space)

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()

    Plot.subspace!(ax, 1:2, sub_space, fa = 0.3)
    Plot.box!(ax, 1:2, [1.0, 0.0], [8.0, 10.0])
    Plot.box!(ax, 1:2, [5.0, 5.0], [10000.0, 10000.0])
    Plot.box!(ax, 1:2, [1.0, 1.0], [2.0, 2.0])
    Plot.box!(ax, 1:2, [0.0, 0.0], [5.0, 5.0])

    ax.set_xlim([-2.0, 14.0])
    ax.set_ylim([-2.0, 14.0])
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
