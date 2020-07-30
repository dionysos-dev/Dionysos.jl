include("../abstraction.jl")

module TestMain

import Main.Abstraction
using PyPlot
AB = Main.Abstraction

function DoTest()
sleep(0.1) # used for good printing
println("Started test")

lb = (0.0, 0.0, -pi - 0.4)
ub = (4.0, 10.0, pi + 0.4)
x0 = (0.0, 0.0, 0.0)
h = (0.2, 0.2, 0.2)
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace_by_box!(X_grid, lb, ub, AB.OUTER)
AB.remove_from_gridspace_by_box!(X_grid, (1.0, 0.0, -100.0), (1.2, 9.0, 100.0), AB.OUTER)
AB.remove_from_gridspace_by_box!(X_grid, (2.2, 0.0, -100.0), (2.4, 5.0, 100.0), AB.OUTER)
AB.remove_from_gridspace_by_box!(X_grid, (2.2, 6.0, -100.0), (2.4, 10.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (3.4, 0.0, -100.0), (3.6, 9.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (4.6, 1.0, -100.0), (4.8, 10.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (5.8, 0.0, -100.0), (6.0, 6.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (5.8, 7.0, -100.0), (6.0, 10.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (7.0, 1.0, -100.0), (7.2, 10.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (8.2, 0.0, -100.0), (8.4, 8.5, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (8.4, 8.2, -100.0), (9.3, 8.6, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (9.3, 7.0, -100.0), (10.0, 7.4, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (8.4, 5.8, -100.0), (9.3, 6.2, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (9.3, 4.6, -100.0), (10.0, 5.0, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (8.4, 3.4, -100.0), (9.3, 3.8, 100.0), AB.OUTER)
# AB.remove_from_gridspace_by_box!(X_grid, (9.3, 2.2, -100.0), (10.0, 2.6, 100.0), AB.OUTER)

X_full = AB.NewSubSet(X_grid)
AB.add_to_subset_all!(X_full)

X_init = AB.NewSubSet(X_grid)
AB.add_to_subset_by_box!(X_init, (0.4, 0.4, 0.0), (0.4, 0.4, 0.0), AB.OUTER)
X_reach = AB.NewSubSet(X_grid)
AB.add_to_subset_by_box!(X_reach, (3.0, 0.5, -100.0), (3.6, 0.8, 100.0), AB.OUTER)

fig = PyPlot.figure()
ax = fig.gca(aspect = "equal")
ax.set_xlim((-0.2, 10.2))
ax.set_ylim((-0.2, 10.2))

AB.plot_subset!(ax, 1:2, X_full, fa = 0.0)
AB.plot_subset!(ax, 1:2, X_init, fc = "green")
AB.plot_subset!(ax, 1:2, X_reach, fc = "yellow")

lb = (-1.0, -1.0)
ub = (1.0, 1.0)
u0 = (0.0, 0.0)
h = (0.3, 0.3)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace_by_box!(U_grid, lb, ub, AB.OUTER)

tstep = 0.3
n_sys = 3
n_bound = 3
function F_sys(x, u)
      alpha = atan(tan(u[2])/2)
      return (
            u[1]*cos(alpha + x[3])/cos(alpha),
            u[1]*sin(alpha + x[3])/cos(alpha),
            u[1]*tan(u[2]))
end
function L_bound(r, u)
      alpha = atan(tan(u[2])/2)
      return (u[1]/cos(alpha)*r[3], u[1]/cos(alpha)*r[3], 0.0)
end
sys_noise = (0.0, 0.0, 0.0)
meas_noise = (0.0, 0.0, 0.0)

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
sym_model_sys = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
@time AB.set_symmodel_from_controlsystem!(sym_model_sys, cont_sys)

sym_model_contr = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
@time AB.set_controller_reach!(sym_model_contr, sym_model_sys, X_init, X_reach)

x0 = (0.4, 0.4, 0.0)
AB.plot_trajectory_closed_loop!(ax, 1:2, cont_sys, sym_model_contr, x0, 100)

x_ref = iterate(AB.enumerate_subset_ref(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(X_grid, x_ref)
X_simple = AB.NewSubSet(X_grid)
XUY_simple_ = Any[]

# for i = 1:6
#     # global x_ref, X_grid, U_grid, XUY_simple_
#     Xs = AB.NewSubSet(X_grid)
#     Ys = AB.NewSubSet(X_grid)
#     Us = AB.NewSubSet(U_grid)
#     AB.add_inputs_images_by_xref!(Us, Xs, sym_model_contr, x_ref)
#     for u_ref in AB.enumerate_subset_ref(Us)
#         AB.add_images_by_xref_uref!(Ys, sym_model_sys, x_ref, u_ref)
#     end
#     push!(XUY_simple_, (Xs, Us, Ys))
#     if ~AB.is_subset_empty(Ys)
#         x_ref = iterate(AB.enumerate_subset_ref(Ys))[1]
#     end
# end
#
# for (Xs, Us, Ys) in XUY_simple_
#     AB.plot_subset!(ax, 1:2, Xs, fc = "green")
#     AB.plot_subset!(ax, 1:2, Ys, fc = "blue")
#     AB.plot_cell_image!(ax, 1:2, Xs, Us, cont_sys)
#     AB.plot_cell_approx!(ax, 1:2, Xs, Us, cont_sys)
# end

sleep(0.1) # used for good printing
println("End test")
end

DoTest()

end  # module TestMain
