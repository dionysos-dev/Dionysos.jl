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
AB.remove_from_gridspace!(X_grid, AB.HyperRectangle((-1.0, -2.0), (-1.1, 4.0)), AB.OUTER)
X_full = AB.NewSubSet(X_grid)
AB.add_to_subset_all!(X_full)

lb = (-2.0,)
ub = (2.0,)
u0 = (0.0,)
h = (1.0,)
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

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

X_init = AB.NewSubSet(X_grid)
AB.add_to_subset!(X_init, AB.HyperRectangle((-3.0, -3.0), (-2.9, -2.9)), AB.OUTER)
X_target = AB.NewSubSet(X_grid)
AB.add_to_subset!(X_target, AB.HyperRectangle((0.0, 0.0), (4.0, 4.0)), AB.OUTER)

sym_model_contr = AB.NewSymbolicModelHash(X_grid, U_grid, X_grid)
AB.set_controller_reach!(sym_model_contr, sym_model_sys, X_init, X_target)
@test AB.get_symmodel_size(sym_model_contr) == 836

x_ref = iterate(AB.enumerate_subset_ref(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(X_grid, x_ref)
X_simple = AB.NewSubSet(X_grid)
XUY_simple_ = Any[]

for i = 1:6
    # global x_ref, X_grid, U_grid, XUY_simple_
    Xs = AB.NewSubSet(X_grid)
    Ys = AB.NewSubSet(X_grid)
    Us = AB.NewSubSet(U_grid)
    xref_coll = AB.get_gridspace_reftype(X_grid)[]
    uref_coll = AB.get_gridspace_reftype(U_grid)[]
    yref_coll = AB.get_gridspace_reftype(X_grid)[]
    AB.add_inputs_images_by_xref!(uref_coll, xref_coll, sym_model_contr, x_ref)
    for u_ref in uref_coll
        AB.add_images_by_xref_uref!(yref_coll, sym_model_sys, x_ref, u_ref)
    end
    AB.add_to_subset_by_ref_coll!(Xs, xref_coll)
    AB.add_to_subset_by_ref_coll!(Us, uref_coll)
    AB.add_to_subset_by_ref_coll!(Ys, yref_coll)
    push!(XUY_simple_, (Xs, Us, Ys))
    if ~AB.is_subset_empty(Ys)
        x_ref = iterate(AB.enumerate_subset_ref(Ys))[1]
    end
end

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-5.3, 5.3))
    ax.set_ylim((-5.3, 5.3))
    Plot.subset!(ax, 1:2, X_full, fa = 0.1)
    Plot.subset!(ax, 1:2, X_init)
    Plot.subset!(ax, 1:2, X_target)
    for (Xs, Us, Ys) in XUY_simple_
        Plot.subset!(ax, 1:2, Xs, fc = "green")
        Plot.subset!(ax, 1:2, Ys, fc = "blue")
        Plot.cell_image!(ax, 1:2, Xs, Us, cont_sys)
        Plot.cell_approx!(ax, 1:2, Xs, Us, cont_sys)
    end
    Plot.trajectory_closed_loop!(ax, 1:2, cont_sys, sym_model_contr, x0, 10)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
