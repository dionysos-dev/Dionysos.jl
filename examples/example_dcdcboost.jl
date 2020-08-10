include("../src/abstraction.jl")
include("../src/plotting.jl")

module DcDcBoost

using Main.Abstraction
using StaticArrays
using PyPlot
using Main.Plot
AB = Main.Abstraction

function dcdc_boost(; nsteps = nothing,
    vs = 1.0, rL = 0.05, xL = 3.0, rC = 0.005, xC = 70.0, r0 = 1.0,
    approx_mode = "nothing")
    #---------------------------------------------------------------------------
    lb = SVector(1.15, 5.45)
    ub = SVector(1.55, 5.85)
    x0 = SVector(0.0, 0.0)
    h = SVector(2.0/4.0e3, 2.0/4.0e3)
    Xgrid = AB.NewGridSpaceList(x0, h)
    AB.add_set!(Xgrid, AB.HyperRectangle(lb, ub), AB.OUTER)
    Xfull = AB.NewSubSet(Xgrid)
    AB.add_all!(Xfull)

    lb = SVector(1)
    ub = SVector(2)
    u0 = SVector(1)
    h = SVector(1)
    Ugrid = AB.NewGridSpaceList(u0, h)
    AB.add_set!(Ugrid, AB.HyperRectangle(lb, ub), AB.OUTER)
    display(Ugrid)

    symmodel = AB.NewSymbolicModelListList(Xgrid, Ugrid)

    Xinit = AB.NewSubSet(Xgrid)
    AB.add_all!(Xinit)
    Xsafe = AB.NewSubSet(Xgrid)
    AB.add_all!(Xsafe)
    initlist = Int[]
    for pos in AB.enum_pos(Xinit)
        push!(initlist, AB.get_state_by_xpos(symmodel, pos))
    end
    safelist = Int[]
    for pos in AB.enum_pos(Xsafe)
        push!(safelist, AB.get_state_by_xpos(symmodel, pos))
    end

    fig = PyPlot.figure()
    ax = fig.gca(aspect = "equal")
    ax.set_xlim((1.15, 1.55))
    ax.set_ylim((5.45, 5.85))

    println("Start plotting")

    # Plot.subset!(ax, 1:2, Xfull, fa = 0.0, ew = 0.5)
    # Plot.subset!(ax, 1:2, Xinit, fc = "green")
    # Plot.subset!(ax, 1:2, Xsafe, fc = "yellow")

    b = SVector(vs/xL, 0.0)
    A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/(xC*(r0+rC)))
    A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, (r0/(r0+rC))/xC,
        -(r0/(r0+rC))/xL, -1.0/(xC*(r0+rC)))
    A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, (r0/(r0+rC))/xC,
        (r0/(r0+rC))/xL, -1.0/(xC*(r0+rC)))
    F_sys = let b = b, A1 = A1, A2 = A2
        function (x, u)
            if u[1] == 1
                return A1*x + b
            elseif u[1] == 2
                return A2*x + b
            end
        end
    end
    L_growthbound = let A1 = A1, A2_abs = A2_abs
        function (u)
            if u[1] == 1
                return A1
            elseif u[1] == 2
                return A2_abs
            end
        end
    end
    DF_sys = let A1 = A1, A2 = A2
        function (x, u)
            if u[1] == 1
                return A1
            elseif u[1] == 2
                return A2
            end
        end
    end
    tstep = 0.05
    nsys = 5
    ngrowthbound = 5
    sysnoise = SVector(0.0, 0.0)
    measnoise = SVector(0.0, 0.0)

    symmodel = AB.NewSymbolicModelListList(Xgrid, Ugrid)

    if approx_mode == "growth"
        contsys = AB.NewControlSystemGrowthRK4(
            tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
    elseif approx_mode == "linearized"
        contsys = AB.NewControlSystemLinearizedRK4(
            tstep, F_sys, DF_sys, bound_DF, bound_DDF, measnoise, nsys)
    end

    @time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

    contr = AB.NewControllerList()
    @time AB.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)

    display(contr)

    x0 = SVector(1.2, 5.6)
    Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, nsteps)
end

end  # module PathPlanning
