include("../src/abstraction.jl")
include("../src/plotting.jl")

module PathPlanning

using Main.Abstraction
using PyPlot
using Main.Plot
AB = Main.Abstraction
using ProfileView

function path_planning(frame_length; nsteps = nothing,
    X1_lb = [1.0, 2.2,  2.2, 3.4,  4.6, 5.8,  5.8,  7.0, 8.2, 8.4,  9.3, 8.4,  9.3, 8.4,  9.3],
    X1_ub = [1.2, 2.4,  2.4, 3.6,  4.8, 6.0,  6.0,  7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0],
    X2_lb = [0.0, 0.0,  6.0, 0.0,  1.0, 0.0,  7.0,  1.0, 0.0, 8.2,  7.0, 5.8,  4.6, 3.4,  2.2],
    X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6,  7.4, 6.2,  5.0, 3.8,  2.6],
    h = (0.2, 0.2, 0.2))
    #---------------------------------------------------------------------------
    frame = AB.HyperRectangle((0.0, 0.0, -pi - 0.4), (frame_length, 10.0, pi + 0.4))
    init = AB.HyperRectangle((0.4, 0.4, 0.0), (0.4, 0.4, 0.0))
    target = AB.HyperRectangle((frame_length - 1.0, 0.5, -100.0), (frame_length - 0.4, 0.8, 100.0))
    x0 = (0.0, 0.0, 0.0)

    Xgrid = AB.NewGridSpaceList(x0, h)
    AB.add_set!(Xgrid, frame, AB.OUTER)
    for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
        box = AB.HyperRectangle((x1lb, x2lb, frame.lb[3]), (x1ub, x2ub, frame.ub[3]))
        if box ⊆ frame && isempty(box ∩ init) && isempty(box ∩ target)
            AB.remove_set!(Xgrid, box, AB.OUTER)
        end
    end
    Xfull = AB.NewSubSet(Xgrid)
    AB.add_all!(Xfull)

    lb = (-1.0, -1.0)
    ub = (1.0, 1.0)
    u0 = (0.0, 0.0)
    h = (0.3, 0.3)
    Ugrid = AB.NewGridSpaceList(u0, h)
    AB.add_set!(Ugrid, AB.HyperRectangle(lb, ub), AB.OUTER)

    symmodel = AB.NewSymbolicModelListList(Xgrid, Ugrid)

    Xinit = AB.NewSubSet(Xgrid)
    AB.add_set!(Xinit, init, AB.OUTER)
    initlist = Int[]
    for pos in AB.enum_pos(Xinit)
        push!(initlist, AB.get_state_by_xpos(symmodel, pos))
    end
    Xtarget = AB.NewSubSet(Xgrid)
    AB.add_set!(Xtarget, target, AB.OUTER)
    targetlist = Int[]
    for pos in AB.enum_pos(Xtarget)
        push!(targetlist, AB.get_state_by_xpos(symmodel, pos))
    end

    fig = PyPlot.figure()
    ax = fig.gca(aspect = "equal")
    ax.set_xlim([-0.2, frame_length + 0.2])
    ax.set_ylim([-0.2, 10.2])

    Plot.subset!(ax, 1:2, Xfull, fa = 0.0)
    Plot.subset!(ax, 1:2, Xinit, fc = "green")
    Plot.subset!(ax, 1:2, Xtarget, fc = "yellow")

    nsteps === nothing && return

    tstep = 0.3
    nsys = 5
    nbound = 5
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
    sysnoise = (0.0, 0.0, 0.0)
    measnoise = (0.0, 0.0, 0.0)

    contsys = AB.NewControlSystemRK4(
        tstep, F_sys, L_bound, sysnoise, measnoise, nsys, nbound)

    @time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

    contr = AB.NewControllerList()
    npoststable = [0 for i = 1:symmodel.autom.nstates, j = 1:symmodel.autom.nsymbols]
    display("npoststable created")
    @time AB.compute_controller_reach!(npoststable, contr, symmodel.autom, initlist, targetlist)

    x0 = (0.4, 0.4, 0.0)
    Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, nsteps)
end

end  # module PathPlanning
