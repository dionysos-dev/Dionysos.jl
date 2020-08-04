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
Xgrid = AB.NewGridSpaceList(x0, h)
AB.add_set!(Xgrid, AB.HyperRectangle(lb, ub), AB.OUTER)
Xfull = AB.NewSubSet(Xgrid)
AB.add_all!(Xfull)

lb = (-4.0,)
ub = (4.0,)
u0 = (0.0,)
h = (0.5,)
Ugrid = AB.NewGridSpaceList(u0, h)
AB.add_set!(Ugrid, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 1.0
nsys = 3
nbound = 3
F_sys(x, u) = (u[1], -x[2] + u[1])
sysnoise = (1.0, 1.0).*0.001
measnoise = (1.0, 1.0).*0.001
L_bound(r, u) = (-0*r[1], -r[2])

contsys = AB.NewControlSystemRK4(tstep, F_sys, L_bound, sysnoise, measnoise, nsys, nbound)
symmodel_sys = AB.NewSymbolicModelListList(Xgrid, Ugrid, Xgrid)
AB.compute_symmodel_from_controlsystem!(symmodel_sys, contsys)

X_init = AB.NewSubSet(Xgrid)
AB.add_set!(X_init, AB.HyperRectangle((-3.0, -3.0), (-2.9, -2.9)), AB.OUTER)

X_safe = AB.NewSubSet(Xgrid)
AB.add_all!(X_safe)
symmodel_contr = AB.NewSymbolicModelListList(Xgrid, Ugrid, Xgrid)
AB.set_controller_safe!(symmodel_contr, symmodel_sys, X_init, X_safe)

correct = true
for x_ref in AB.enum_pos(Xfull)
    if !correct
        break
    end
    uref_coll1 = AB.getgridspace_reftype(Ugrid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll1, symmodel_contr, x_ref, Xfull)
    uref_coll2 = AB.getgridspace_reftype(Ugrid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll2, symmodel_sys, x_ref, Xfull)
    correct = Set(uref_coll1) == Set(uref_coll2)
end
@test correct

AB.remove_fromsubset!(X_safe, AB.HyperRectangle((-1.0, -2.0), (-1.1, 4.0)), AB.OUTER)
symmodel_contr = AB.NewSymbolicModelListList(Xgrid, Ugrid, Xgrid)
AB.set_controller_safe!(symmodel_contr, symmodel_sys, X_init, X_safe)

X_inv = AB.NewSubSet(Xgrid)
Y_safe = AB.NewSubSet(Xgrid)
U_safe = AB.NewSubSet(Ugrid)

for x_ref in AB.enum_pos(X_safe)
    uref_coll = AB.getgridspace_reftype(Ugrid)[]
    AB.add_inputs_by_xref_ysub!(uref_coll, symmodel_contr, x_ref, Xfull)
    AB.add_pos_coll!(U_safe, uref_coll)
    yref_coll = AB.getgridspace_reftype(Xgrid)[]
    for u_ref in uref_coll
        AB.add_images_by_xref_uref!(yref_coll, symmodel_sys, x_ref, u_ref)
    end
    AB.add_pos_coll!(Y_safe, yref_coll)
    if !isempty(uref_coll)
        AB.add_set_by_new_ref!(X_inv, x_ref)
    end
end

x_ref = iterate(AB.enum_pos(X_init))[1]
display(x_ref)
x0 = AB.get_coords_by_ref(Xgrid, x_ref)

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-5.5, 5.5))
    ax.set_ylim((-5.3, 5.3))
    Plot.subset!(ax, 1:2, Xfull, fa = 0.0)
    Plot.subset!(ax, 1:2, X_init)
    Plot.subset!(ax, 1:2, X_safe, fa = 0.1)
    Plot.subset!(ax, 1:2, X_inv, fa = 0.1, fc = "yellow")
    Plot.subset!(ax, 1:2, Y_safe; fc = "blue")
    Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel_contr, x0, 100)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
