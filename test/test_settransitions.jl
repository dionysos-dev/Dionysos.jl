include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays

import Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "Transitions" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace!(X_grid, AB.HyperRectangle(lb, ub), AB.OUTER)
X_full = AB.NewSubSpace(X_grid)
AB.add_to_subspace_all!(X_full)

lb = SVector(-1.0)
ub = SVector(1.0)
u0 = SVector(0.0)
h = SVector(0.5)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 5.0
n_sys = 3
n_bound = 3
# F_sys(x, u) = [1.0-cos(x[2]), -x[1] + u[1]]
# L_bound(u) = [0.0 1.0; 1.0 0.0]
F_sys(x, u) = SVector(u[1], -cos(x[1]))
L_bound(u) = @SMatrix[0.0 0.0; 1.0 0.0]
sys_noise = SVector(1.0, 1.0)*0.1
meas_noise = SVector(1.0, 1.0)*0.01*0

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
trans_map = trans_map = AB.NewTransitionMapHash(X_grid, U_grid, X_grid)
AB.set_transitions_from_controlsystem!(trans_map, cont_sys)
@test length(trans_map.elems) == 1145

x_pos = SVector(1, 2)
u_pos = SVector(1)
x = AB.get_coords_by_pos(X_grid, x_pos)
u = AB.get_coords_by_pos(U_grid, u_pos)
x_ref = AB.get_ref_by_pos(X_grid, x_pos)
u_ref = AB.get_ref_by_pos(U_grid, u_pos)
y_ref_coll = AB.get_transition_image(trans_map, x_ref, u_ref)[1][2]
@test length(y_ref_coll) == 18
display(y_ref_coll)

X_simple = AB.NewSubSpace(X_grid)
AB.add_to_subspace_by_pos!(X_simple, x_pos)
U_simple = AB.NewSubSpace(U_grid)
AB.add_to_subspace_by_pos!(U_simple, u_pos)
Y_simple = AB.NewSubSpace(X_grid)
for y_ref in y_ref_coll
    AB.add_to_subspace_by_ref!(Y_simple, y_ref)
end

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()

    Plot.subspace!(ax, 1:2, X_full, fa = 0.1)
    Plot.subspace!(ax, 1:2, X_simple)
    Plot.subspace!(ax, 1:2, Y_simple)
    Plot.trajectory_open_loop!(ax, 1:2, cont_sys, x, u, 50)
    Plot.cell_image!(ax, 1:2, X_simple, U_simple, cont_sys)
    Plot.cell_approx!(ax, 1:2, X_simple, U_simple, cont_sys)

    ax.set_xlim([-1.0, 11.0])
    ax.set_ylim([-2.0, 14.0])
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
