include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))
include(joinpath(@__DIR__, "../../src/Abstraction/plotting.jl"))

module PathPlanning

using Main.Abstraction
using Main.Plot
using StaticArrays
using PyPlot
ABS = Main.Abstraction

include("pathplanning_system.jl")

function path_planning(x1_max; nstep = nothing,
    X1_lb = [1.0, 2.2,  2.2, 3.4,  4.6, 5.8,  5.8,  7.0, 8.2, 8.4,  9.3, 8.4,  9.3, 8.4,  9.3],
    X1_ub = [1.2, 2.4,  2.4, 3.6,  4.8, 6.0,  6.0,  7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0],
    X2_lb = [0.0, 0.0,  6.0, 0.0,  1.0, 0.0,  7.0,  1.0, 0.0, 8.2,  7.0, 5.8,  4.6, 3.4,  2.2],
    X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6,  7.4, 6.2,  5.0, 3.8,  2.6],
    xstep = SVector(0.2, 0.2, 0.2), approx_mode = "nothing")
    #---------------------------------------------------------------------------
    xorig = SVector(0.0, 0.0, 0.0)
    xh = xstep
    uorig = SVector(0.0, 0.0)
    uh = SVector(0.3, 0.3)
    mng = ABS.ListGridManager(xorig, xh, uorig, uh)
    XDom = ABS.AddXDomain!(mng)
    rectFullX = ABS.HyperRectangle(SVector(0.0, 0.0, -pi-0.4), SVector(x1_max, 10.0, pi+0.4))
    rectInit = ABS.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
    rectTarget = ABS.HyperRectangle(SVector(x1_max-1.0, 0.5, -Inf), SVector(x1_max-0.4, 0.8, Inf))
    ABS.add_cells!(mng, XDom, rectFullX, ABS.OUTER)
    for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
        rect = ABS.HyperRectangle(SVector(x1lb, x2lb, -Inf), SVector(x1ub, x2ub, Inf))
        if isempty(rect ∩ rectInit) && isempty(rect ∩ rectTarget)
            ABS.remove_cells!(mng, XDom, rect, ABS.OUTER)
        end
    end
    ABS.create_symbols!(mng, XDom)
    stateset = ABS.AddStateSet!(mng)
    ABS.add_symbols!(mng, stateset, XDom)

    println("XDom and states created")

    UDom = ABS.AddUDomain!(mng)
    rectFullU = ABS.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))
    ABS.add_cells!(mng, UDom, rectFullU, ABS.OUTER)
    ABS.create_symbols!(mng, UDom)
    labelset = ABS.AddLabelSet!(mng)
    ABS.add_symbols!(mng, labelset, UDom)

    println("UDom and labels created")

    XInitDom = ABS.AddXDomain!(mng)
    ABS.add_cells!(mng, XInitDom, XDom, rectInit, ABS.OUTER)
    initset = ABS.AddStateSet!(mng)
    ABS.add_symbols!(mng, initset, XInitDom)

    XTargetDom = ABS.AddXDomain!(mng)
    ABS.add_cells!(mng, XTargetDom, XDom, rectTarget, ABS.OUTER)
    targetset = ABS.AddStateSet!(mng)
    ABS.add_symbols!(mng, targetset, XTargetDom)

    println("Init and target sets created")

    if approx_mode == "growth"
        contsys = ABS.ControlSystemGrowthRK4(
            tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
    elseif approx_mode == "linearized"
        contsys = ABS.ControlSystemLinearizedRK4(
            tstep, F_sys, DF_sys, bound_DF, bound_DDF, measnoise, nsys)
    end

    println("Control system created")

    autom = ABS.AddAutomaton!(mng)
    @time ABS.add_transitions!(mng, autom, stateset, labelset, stateset, contsys)
    empty!(mng, autom)
    @time ABS.add_transitions!(mng, autom, stateset, labelset, stateset, contsys)

    println("Automaton created")

    spec = ABS.ReachSpec(mng, initset, targetset)
    contr = ABS.AddController!(mng)
    @time ABS.add_controls!(mng, contr, autom, spec)
    empty!(mng, contr)
    @time ABS.add_controls!(mng, contr, autom, spec)

    if nstep !== nothing
        fig = PyPlot.figure()
        ax = fig.gca(aspect = "equal")
        ax.set_xlim((-0.2, x1_max + 0.2))
        ax.set_ylim((-0.2, 10.2))

        Plot.cells!(ax, 1:2, mng, XDom, fa = 0.0)
        Plot.cells!(ax, 1:2, mng, XInitDom, fc = "green")
        Plot.cells!(ax, 1:2, mng, XTargetDom, fc = "yellow")

        x0 = SVector(0.4, 0.4, 0.0)
        Plot.trajectory!(ax, 1:2, mng, contr, contsys, x0, nstep)
    end
end

end  # module PathPlanning
