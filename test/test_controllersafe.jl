include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControllerReach" begin
lb = (-5.0, -5.0)
ub = (5.0, 5.0)
x0 = (0.0, 0.0)
h = (0.47, 0.23)
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace!(X_grid, AB.HyperRectangle(lb, ub), AB.OUTER)
X_full = AB.NewSubSet(X_grid)
AB.add_to_subset_all!(X_full)

lb = (-4.0,)
ub = (4.0,)
u0 = (0.0,)
h = (0.5,)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 1.0
n_sys = 3
n_bound = 3
F_sys(x, u) = (u[1], -x[2] + u[1])
sys_noise = (1.0, 1.0).*0.001
meas_noise = (1.0, 1.0).*0.001
L_bound(r, u) = (-0*r[1], -r[2])

cont_sys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
sym_model_sys = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_symmodel_from_controlsystem!(sym_model_sys, cont_sys)

X_init = AB.NewSubSet(X_grid)
AB.add_to_subset!(X_init, AB.HyperRectangle((-3.0, -3.0), (-2.9, -2.9)), AB.OUTER)

X_safe = AB.NewSubSet(X_grid)
AB.add_to_subset_all!(X_safe)
sym_model_contr = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_controller_safe!(sym_model_contr, sym_model_sys, X_init, X_safe)

correct = true
for x_ref in AB.enumerate_subset_ref(X_full)
    if !correct
        break
    end
    uref_coll1 = AB.get_gridspace_reftype(U_grid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll1, sym_model_contr, x_ref, X_full)
    uref_coll2 = AB.get_gridspace_reftype(U_grid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll2, sym_model_sys, x_ref, X_full)
    correct = Set(uref_coll1) == Set(uref_coll2)
end
@test correct

AB.remove_from_subset!(X_safe, AB.HyperRectangle((-1.0, -2.0), (-1.1, 4.0)), AB.OUTER)
sym_model_contr = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_controller_safe!(sym_model_contr, sym_model_sys, X_init, X_safe)

X_inv = AB.NewSubSet(X_grid)
Y_safe = AB.NewSubSet(X_grid)
U_safe = AB.NewSubSet(U_grid)

for x_ref in AB.enumerate_subset_ref(X_safe)
    uref_coll = AB.get_gridspace_reftype(U_grid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll, sym_model_contr, x_ref, X_full)
    AB.add_to_subset_by_ref_coll!(U_safe, uref_coll)
    yref_coll = AB.get_gridspace_reftype(X_grid)[]
    for u_ref in uref_coll
        AB.add_images_by_xref_uref!(yref_coll, sym_model_sys, x_ref, u_ref)
    end
    AB.add_to_subset_by_ref_coll!(Y_safe, yref_coll)
    if !isempty(uref_coll)
        AB.add_to_subset_by_new_ref!(X_inv, x_ref)
    end
end

x_ref = iterate(AB.enumerate_subset_ref(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(X_grid, x_ref)

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-5.5, 5.5))
    ax.set_ylim((-5.3, 5.3))
    Plot.subset!(ax, 1:2, X_full, fa = 0.0)
    Plot.subset!(ax, 1:2, X_init)
    Plot.subset!(ax, 1:2, X_safe, fa = 0.1)
    Plot.subset!(ax, 1:2, X_inv, fa = 0.1, fc = "yellow")
    Plot.subset!(ax, 1:2, Y_safe; fc = "blue")
    Plot.trajectory_closed_loop!(ax, 1:2, cont_sys, sym_model_contr, x0, 100)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
