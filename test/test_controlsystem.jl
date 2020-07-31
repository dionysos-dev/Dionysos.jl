include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControlSystem" begin
lb = (0.0, 0.0)
ub = (10.0, 11.0)
x0 = (0.0, 0.0)
h = (1.0, 2.0)
X_grid = AB.NewGridSpaceHash(x0, h)
@test AB.get_gridspace_size(X_grid) == 0
AB.add_to_gridspace!(X_grid, AB.HyperRectangle(lb, ub), AB.OUTER)
@test AB.get_gridspace_size(X_grid) == 77

lb = (-1.0,)
ub = (1.0,)
u0 = (0.0,)
h = (0.5,)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)
@test AB.get_gridspace_size(U_grid) == 5

tstep = 1.0
n_sys = 3
n_bound = 3
# F_sys(x, u) = (1.0-cos(x(2)), -x(1) + u(1)
# L_bound(u) = (0.0 1.0; 1.0 0.0)
F_sys(x, u) = (u[1], -cos(x[1]))
L_bound(r, u) = (0.0, r[1])
sys_noise = (1.0, 1.0).*0.01
meas_noise = (1.0, 1.0).*0.01

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)


x_pos = (1, 1)
u_pos = (1,)
x = AB.get_coords_by_pos(X_grid, x_pos)
u = AB.get_coords_by_pos(U_grid, u_pos)
X_simple = AB.NewSubSet(X_grid)
AB.add_to_subset_by_pos!(X_simple, x_pos)
@test AB.get_subset_size(X_simple) == 1
U_simple = AB.NewSubSet(U_grid)
AB.add_to_subset_by_pos!(U_simple, u_pos)
@test AB.get_subset_size(U_simple) == 1

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-10.0, 10.0))
    ax.set_ylim((-10.0, 10.0))
    Plot.subset!(ax, 1:2, X_simple)
    Plot.trajectory_open_loop!(ax, 1:2, cont_sys, x, u, 50)
    Plot.cell_image!(ax, 1:2, X_simple, U_simple, cont_sys)
    Plot.cell_approx!(ax, 1:2, X_simple, U_simple, cont_sys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
