include("../src/abstraction.jl")
include("../src/plotting.jl")

module PathPlanning

using Main.Abstraction
using StaticArrays
using PyPlot
using Main.Plot
AB = Main.Abstraction

function path_planning(x1_max; nstep = nothing,
    X1_lb = [1.0, 2.2,  2.2, 3.4,  4.6, 5.8,  5.8,  7.0, 8.2, 8.4,  9.3, 8.4,  9.3, 8.4,  9.3],
    X1_ub = [1.2, 2.4,  2.4, 3.6,  4.8, 6.0,  6.0,  7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0],
    X2_lb = [0.0, 0.0,  6.0, 0.0,  1.0, 0.0,  7.0,  1.0, 0.0, 8.2,  7.0, 5.8,  4.6, 3.4,  2.2],
    X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6,  7.4, 6.2,  5.0, 3.8,  2.6],
    h = SVector(0.2, 0.2, 0.2), approx_mode = "nothing")
    #---------------------------------------------------------------------------
    _X_ = AB.HyperRectangle(SVector(0.0, 0.0, -pi-0.4), SVector(x1_max, 10.0, pi+0.4))
    _I_ = AB.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
    _T_ = AB.HyperRectangle(SVector(x1_max-1.0, 0.5, -100.0), SVector(x1_max-0.4, 0.8, 100.0))
    x0 = SVector(0.0, 0.0, 0.0)
    Xgrid = AB.GridFree(x0, h)
    Xfull = AB.DomainList(Xgrid)
    AB.add_set!(Xfull, _X_, AB.OUTER)
    for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
        box = AB.HyperRectangle(SVector(x1lb, x2lb, _X_.lb[3]), SVector(x1ub, x2ub, _X_.ub[3]))
        if box ⊆ _X_ && isempty(box ∩ _I_) && isempty(box ∩ _T_)
            AB.remove_set!(Xfull, box, AB.OUTER)
        end
    end

    _U_ = AB.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))
    u0 = SVector(0.0, 0.0)
    h = SVector(0.3, 0.3)
    Ugrid = AB.GridFree(u0, h)
    Ufull = AB.DomainList(Ugrid)
    AB.add_set!(Ufull, _U_, AB.OUTER)

    symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)

    Xinit = AB.DomainList(Xgrid)
    AB.add_subset!(Xinit, Xfull, _I_, AB.OUTER)
    initlist = Int[]
    for pos in AB.enum_pos(Xinit)
        push!(initlist, AB.get_state_by_xpos(symmodel, pos))
    end
    Xtarget = AB.DomainList(Xgrid)
    AB.add_subset!(Xtarget, Xfull, _T_, AB.OUTER)
    targetlist = Int[]
    for pos in AB.enum_pos(Xtarget)
        push!(targetlist, AB.get_state_by_xpos(symmodel, pos))
    end

    #---------------------------------------------------------------------------
    tstep = 0.3
    nsys = 5
    ngrowthbound = 5
    function F_sys(x, u)
        α = atan(tan(u[2])/2)
        return SVector{3}(
            u[1]*cos(α + x[3])/cos(α),
            u[1]*sin(α + x[3])/cos(α),
            u[1]*tan(u[2]))
    end
    # There was a mistake below: forgot the abs !!!
    # Now we have 23784452 transitions for pathplanning-hard
    function L_growthbound(u)
        β = abs(u[1]/cos(atan(tan(u[2])/2)))
        # Both have the same speed
        # return @SMatrix [
        #     0.0 0.0 β;
        #     0.0 0.0 β;
        #     0.0 0.0 0.0]
        return SMatrix{3,3}(
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            β, β, 0.0)
    end
    function DF_sys(x, u)
        α = atan(tan(u[2])/2)
        β = u[1]/cos(α)
        return SMatrix{3,3}(
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                -β*sin(α + x[3]), β*cos(α + x[3]), 0.0)
    end
    bound_DF(u) = abs(u[1]/cos(atan(tan(u[2])/2)))
    # DDF_1 = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -β*cos(α + x[3])]
    # DDF_2 = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -β*sin(α + x[3])]
    # DDF_3 = 0.0
    bound_DDF(u) = abs(u[1]/cos(atan(tan(u[2])/2)))
    sysnoise = SVector(0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0)

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
    # Uncomment for fair eval of perfs
    empty!(symmodel.autom)
    @time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

    # Uncomment to compare new and OLD versions
    # empty!(symmodel.autom)
    # @time AB.compute_symmodel_from_controlsystem_OLD!(symmodel, contsys)
    # empty!(symmodel.autom)
    # @time AB.compute_symmodel_from_controlsystem_OLD!(symmodel, contsys)

    contr = AB.NewControllerList()

    empty!(contr)
    @time AB.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)
    # Uncomment for fair eval of perfs
    empty!(contr)
    @time AB.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)


    if nstep !== nothing
        fig = PyPlot.figure()
        ax = fig.gca(aspect = "equal")
        ax.set_xlim((-0.2, x1_max + 0.2))
        ax.set_ylim((-0.2, 10.2))

        Plot.domain!(ax, 1:2, Xfull, fa = 0.0)
        Plot.domain!(ax, 1:2, Xinit, fc = "green")
        Plot.domain!(ax, 1:2, Xtarget, fc = "yellow")

        x0 = SVector(0.4, 0.4, 0.0)
        Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, nstep)
    end
end

end  # module PathPlanning
