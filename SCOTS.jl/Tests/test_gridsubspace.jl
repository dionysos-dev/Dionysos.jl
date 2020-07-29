include("../abstraction.jl")

module TestMain

using Main.Abstraction
using PyPlot
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

orig = [0.0, 0.0]
h = [1.0, 2.0]
grid_space = AB.NewGridSpaceHash(orig, h)
display(grid_space)

AB.add_to_gridspace_by_coords!(grid_space, [1.2, 3.5])
AB.add_to_gridspace_by_coords!(grid_space, [-0.5002, -1.0])
display(grid_space)

AB.add_to_gridspace_by_box!(grid_space, [1.0, 0.0], [11.0, 10.0], AB.OUTER)
println(grid_space.elems)

AB.remove_from_gridspace_by_coords!(grid_space, [2.0, 2.0])
AB.remove_from_gridspace_by_box!(grid_space, [5.0, 5.0], [10000.0, 10000.0], AB.INNER)
println(grid_space.elems)

pos_coll = AB.enumerate_gridspace_pos(grid_space)
display(pos_coll)

sub_space = AB.NewSubSpace(grid_space)
AB.add_to_subspace_all!(sub_space)
AB.remove_from_subspace_by_box!(sub_space, [1.0, 1.0], [2.0, 2.0], AB.OUTER)
AB.add_to_subspace_by_box!(sub_space, [0.0, 0.0], [5.0, 5.0], AB.INNER)
display(sub_space)

fig = PyPlot.figure()
ax = fig.gca()

AB.plot_subspace!(ax, 1:2, sub_space, fa = 0.3)
AB.plot_box!(ax, 1:2, [1.0, 0.0], [8.0, 10.0])
AB.plot_box!(ax, 1:2, [5.0, 5.0], [10000.0, 10000.0])
AB.plot_box!(ax, 1:2, [1.0, 1.0], [2.0, 2.0])
AB.plot_box!(ax, 1:2, [0.0, 0.0], [5.0, 5.0])

ax.set_xlim([-2.0, 14.0])
ax.set_ylim([-2.0, 14.0])

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
