include("../abstraction.jl")

module TestMain

import Main.Abstraction
using PyPlot
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

lb = (-5.0, -5.0)
ub = (5.0, 5.0)
x0 = (0.0, 0.0)
h = (0.47, 0.23)
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace_by_box!(X_grid, lb, ub, AB.INNER)
AB.remove_from_gridspace_by_box!(X_grid, (-1.0, -2.0), (-1.1, 4.0), AB.OUTER)
X_full = AB.NewSubSpace(X_grid)
AB.add_to_subspace_all!(X_full)

lb = (-2.0,)
ub = (2.0,)
u0 = (0.0,)
h = (1.0,)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace_by_box!(U_grid, lb, ub, AB.OUTER)

tstep = 1.0
n_sys = 3
n_bound = 3
F_sys(x, u) = (1.0, u[1])
sys_noise = (1.0, 1.0).*0.001
meas_noise = (1.0, 1.0).*0.001
L_bound(r, u) = (0.0, 0.0)

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
sym_model_sys = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_symmodel_from_controlsystem!(sym_model_sys, cont_sys)

X_init = AB.NewSubSpace(X_grid)
AB.add_to_subspace_by_box!(X_init, (-3.0, -3.0), (-2.9, -2.9), AB.OUTER)
X_reach = AB.NewSubSpace(X_grid)
AB.add_to_subspace_by_box!(X_reach, (0.0, 0.0), (4.0, 4.0), AB.OUTER)

fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-5.3, 5.3))
ax.set_ylim((-5.3, 5.3))

sym_model_contr = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_controller_reach!(sym_model_contr, sym_model_sys, X_init, X_reach)

x_ref = iterate(AB.enumerate_subspace_ref(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(X_grid, x_ref)
X_simple = AB.NewSubSpace(X_grid)
XUY_simple_ = Any[]

for i = 1:6
    global x_ref, X_grid, U_grid, XUY_simple_
    Xs = AB.NewSubSpace(X_grid)
    Ys = AB.NewSubSpace(X_grid)
    Us = AB.NewSubSpace(U_grid)
    AB.add_inputs_images_by_xref!(Us, Xs, sym_model_contr, x_ref)
    for u_ref in AB.enumerate_subspace_ref(Us)
        AB.add_images_by_xref_uref!(Ys, sym_model_sys, x_ref, u_ref)
    end
    push!(XUY_simple_, (Xs, Us, Ys))
    if ~AB.is_subspace_empty(Ys)
        x_ref = iterate(AB.enumerate_subspace_ref(Ys))[1]
    end
end

AB.plot_subspace!(ax, 1:2, X_full, fa = 0.1)
AB.plot_subspace!(ax, 1:2, X_init)
AB.plot_subspace!(ax, 1:2, X_reach)

for (Xs, Us, Ys) in XUY_simple_
    AB.plot_subspace!(ax, 1:2, Xs, fc = "green")
    AB.plot_subspace!(ax, 1:2, Ys, fc = "blue")
    AB.plot_cell_image!(ax, 1:2, Xs, Us, cont_sys)
    AB.plot_cell_approx!(ax, 1:2, Xs, Us, cont_sys)
end

AB.plot_trajectory_closed_loop!(ax, 1:2, cont_sys, sym_model_contr, x0, 10)

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
