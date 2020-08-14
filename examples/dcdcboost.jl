include("../src/abstraction.jl")
include("../src/plotting.jl")

module DcdcBoost

using Main.Abstraction
using StaticArrays
using PyPlot
using Main.Plot
AB = Main.Abstraction

function dcdc_boost(; nstep = nothing,
    vs = 1.0, rL = 0.05, xL = 3.0, rC = 0.005, xC = 70.0, r0 = 1.0,
    approx_mode = "nothing")
    #---------------------------------------------------------------------------
    _X_ = AB.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85))
    x0 = SVector(0.0, 0.0)
    h = SVector(2.0/4.0e3, 2.0/4.0e3)
    Xgrid = AB.GridFree(x0, h)
    Xfull = AB.DomainList(Xgrid)
    AB.add_set!(Xfull, _X_, AB.INNER)

    _U_ = AB.HyperRectangle(SVector(1), SVector(2))
    u0 = SVector(1)
    h = SVector(1)
    Ugrid = AB.GridFree(u0, h)
    Ufull = AB.DomainList(Ugrid)
    AB.add_set!(Ufull, _U_, AB.OUTER)

    symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)

    Xinit = AB.DomainList(Xgrid)
    union!(Xinit, Xfull)
    Xsafe = AB.DomainList(Xgrid)
    union!(Xsafe, Xfull)
    initlist = Int[]
    for pos in AB.enum_pos(Xinit)
        push!(initlist, AB.get_state_by_xpos(symmodel, pos))
    end
    safelist = Int[]
    for pos in AB.enum_pos(Xsafe)
        push!(safelist, AB.get_state_by_xpos(symmodel, pos))
    end

    #---------------------------------------------------------------------------
    b = SVector(vs/xL, 0.0)
    A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC))
    A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
        -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
        r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    F_sys = let b = b, A1 = A1, A2 = A2
        (x, u) -> u[1] == 1 ? A1*x + b : A2*x + b
    end
    L_growthbound = let A1 = A1, A2_abs = A2_abs
        u -> u[1] == 1 ? A1 : A2_abs
    end
    DF_sys = let A1 = A1, A2 = A2
        (x, u) -> u[1] == 1 ? A1 : A2
    end
    bound_DF(u) = 0.0
    bound_DDF(u) = 0.0
    tstep = 0.5
    nsys = 5
    ngrowthbound = 5
    sysnoise = SVector(0.0, 0.0)
    measnoise = SVector(0.0, 0.0)

    if approx_mode == "growth"
        contsys = AB.NewControlSystemGrowthRK4(
            tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
    elseif approx_mode == "linearized"
        contsys = AB.NewControlSystemLinearizedRK4(
            tstep, F_sys, DF_sys, bound_DF, bound_DDF, measnoise, nsys)
    end
    #---------------------------------------------------------------------------

    empty!(symmodel.autom)
    @time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)
    empty!(symmodel.autom)
    @time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

    contr = AB.NewControllerList()

    empty!(contr)
    @time AB.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)
    empty!(contr)
    @time AB.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)

    if nstep !== nothing
        fig = PyPlot.figure()
        ax = fig.gca(aspect = "equal")
        ax.set_xlim((1.15, 1.55))
        ax.set_ylim((5.45, 5.85))

        # Plot.domain!(ax, 1:2, Xfull, fa = 0.0, ew = 0.5)
        # Plot.domain!(ax, 1:2, Xinit, fc = "green")
        # Plot.domain!(ax, 1:2, Xsafe, fc = "yellow")

        # xpos = AB.get_pos_by_coord(Xgrid, SVector(1.2, 5.6))
        # upos = AB.get_pos_by_coord(Ugrid, SVector(2))
        # x = AB.get_coord_by_pos(Xgrid, xpos)
        # u = AB.get_coord_by_pos(Ugrid, upos)
        # source = AB.get_state_by_xpos(symmodel, xpos)
        # symbol = AB.get_symbol_by_upos(symmodel, upos)
        # Xsimple = AB.DomainList(Xgrid)
        # AB.add_pos!(Xsimple, xpos)
        # Usimple = AB.DomainList(Ugrid)
        # AB.add_pos!(Usimple, upos)
        # Ysimple = AB.DomainList(Xgrid)
        # targetlist = Int[]
        # AB.compute_post!(targetlist, symmodel.autom, source, symbol)
        # for target in targetlist
        #     AB.add_pos!(Ysimple, AB.get_xpos_by_state(symmodel, target))
        # end
        #
        # Plot.domain!(ax, 1:2, Xsimple)
        # Plot.domain!(ax, 1:2, Ysimple)
        # Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
        # Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
        # Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)

        x0 = SVector(1.2, 5.6)
        Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, nstep)
    end
end

end  # module PathPlanning
