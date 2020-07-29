include("../src/abstraction.jl")
include("../src/plotting.jl")

module TestMain

import Main.Abstraction
using PyPlot
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

lb = [-5.0, -5.0]
ub = [5.0, 5.0]
x0 = [0.0, 0.0]
h = [0.47, 0.23]
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace!(X_grid, AB.HyperRectangle(lb, ub), AB.INNER)
X_full = AB.NewSubSpace(X_grid)
AB.add_to_subspace_all!(X_full)

lb = [-2.0]
ub = [2.0]
u0 = [0.0]
h = [1.0]
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 1.0
n_sys = 3
n_bound = 3
F_sys(x, u) = [1.0, u[1]]
sys_noise = [1.0, 1.0]*0.001
meas_noise = [1.0, 1.0]*0.001
L_bound(u) = [0.0 0.0; 0.0 0.0]

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
trans_map_sys = AB.NewTransitionMapHash(X_grid, U_grid, X_grid)
AB.set_transitions_from_controlsystem!(trans_map_sys, cont_sys)

X_init = AB.NewSubSpace(X_grid)
AB.add_to_subspace!(X_init, AB.HyperRectangle([-3.0, -3.0], [-2.9, -2.9]), AB.OUTER)
X_reach = AB.NewSubSpace(X_grid)
AB.add_to_subspace!(X_reach, AB.HyperRectangle([0.0, 0.0], [4.0, 4.0]), AB.OUTER)

trans_map_contr = AB.NewTransitionMapHash(X_grid, U_grid, X_grid)
AB.set_controller_reach!(trans_map_contr, trans_map_sys, X_init, X_reach)

x_ref = iterate(AB.enumerate_subspace_ref(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(X_grid, x_ref)
X_simple = AB.NewSubSpace(X_grid)
XUY_simple_ = Any[]

for i = 1:6
    global x_ref, X_grid, U_grid, XUY_simple_
    Xs = AB.NewSubSpace(X_grid)
    AB.add_to_subspace_by_ref!(Xs, x_ref)
    Ys = AB.NewSubSpace(X_grid)
    Us = AB.NewSubSpace(U_grid)
    uy_ref_coll = AB.get_transition_image(trans_map_contr, x_ref)
    if isempty(uy_ref_coll)
        break
    end
    for uy_ref in uy_ref_coll
        AB.add_to_subspace_by_ref!(Us, uy_ref[1])
        for y_ref in uy_ref[2]
            AB.add_to_subspace_by_ref!(Ys, y_ref)
        end
    end
    push!(XUY_simple_, (Xs, Us, Ys))
    x_ref = iterate(uy_ref_coll)[1][2][1]
end

fig = PyPlot.figure()
ax = fig.gca()

Plot.subspace!(ax, 1:2, X_full, fa = 0.1)
Plot.subspace!(ax, 1:2, X_init)
Plot.subspace!(ax, 1:2, X_reach)

for (Xs, Us, Ys) in XUY_simple_
    Plot.subspace!(ax, 1:2, Xs, fc = "green")
    Plot.subspace!(ax, 1:2, Ys, fc = "blue")
    Plot.cell_image!(ax, 1:2, Xs, Us, cont_sys)
    Plot.cell_approx!(ax, 1:2, Xs, Us, cont_sys)
end

Plot.trajectory_closed_loop!(ax, 1:2, cont_sys, trans_map_contr, x0, 10)

ax.set_xlim([-5.3, 5.3])
ax.set_ylim([-5.3, 5.3])

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
