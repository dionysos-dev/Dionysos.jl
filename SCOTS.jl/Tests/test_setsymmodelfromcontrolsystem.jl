include("../abstraction.jl")

module TestMain

import Main.Abstraction
using PyPlot
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

lb = (0.0, 0.0)
ub = (10.0, 11.0)
x0 = (0.0, 0.0)
h = (1.0, 2.0)
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace_by_box!(X_grid, lb, ub, AB.OUTER)
X_full = AB.NewSubSpace(X_grid)
AB.add_to_subspace_all!(X_full)

lb = (-1.0,)
ub = (1.0,)
u0 = (0.0,)
h = (0.5,)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace_by_box!(U_grid, lb, ub, AB.OUTER)

tstep = 5.0
n_sys = 3
n_bound = 3
# F_sys(x, u) = (1.0-cos(x(2)), -x(1) + u(1))
# L_bound(u) = (0.0 1.0; 1.0 0.0)
F_sys(x, u) = (u[1], -cos(x[1]))
L_bound(r, u) = (0.0, r[1])
sys_noise = (1.0, 1.0).*0.0
meas_noise = (1.0, 1.0).*0.0

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
sym_model = sym_model = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_symmodel_from_controlsystem!(sym_model, cont_sys)

x_pos = (1, 2)
u_pos = (1,)
x = AB.get_coords_by_pos(X_grid, x_pos)
u = AB.get_coords_by_pos(U_grid, u_pos)
x_ref = AB.get_ref_by_pos(X_grid, x_pos)
u_ref = AB.get_ref_by_pos(U_grid, u_pos)

X_simple = AB.NewSubSpace(X_grid)
AB.add_to_subspace_by_pos!(X_simple, x_pos)
U_simple = AB.NewSubSpace(U_grid)
AB.add_to_subspace_by_pos!(U_simple, u_pos)
Y_simple = AB.NewSubSpace(X_grid)
yref_vec = AB.add_images_by_xref_uref!(Y_simple, sym_model, x_ref, u_ref)
display(Y_simple)

fig = PyPlot.figure()
ax = fig.gca()

AB.plot_subspace!(ax, 1:2, X_full, fa = 0.1)
AB.plot_subspace!(ax, 1:2, X_simple)
AB.plot_subspace!(ax, 1:2, Y_simple)
AB.plot_trajectory_open_loop!(ax, 1:2, cont_sys, x, u, 50)
AB.plot_cell_image!(ax, 1:2, X_simple, U_simple, cont_sys)
AB.plot_cell_approx!(ax, 1:2, X_simple, U_simple, cont_sys)

ax.set_xlim((-1.0, 11.0))
ax.set_ylim((-2.0, 14.0))

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
