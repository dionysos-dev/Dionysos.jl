include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))
include(joinpath(@__DIR__, "../../src/Abstraction/plotting.jl"))

module TestMain

using Main.Abstraction
ABS = Main.Abstraction
using Main.Plot
using StaticArrays
using PyPlot

xorig = SVector(0.0, 0.0)
xh = SVector(0.47, 0.23)
uorig = SVector(0.0)
uh = SVector(1.0)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

XDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(-5.0, -5.0), SVector(5.0, 5.0))
ABS.add_cells!(mng, XDom, rect, ABS.OUTER)
rect = ABS.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0))
ABS.remove_cells!(mng, XDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, XDom)
stateset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, stateset, XDom)

UDom = ABS.AddUDomain!(mng)
rect = rect = ABS.HyperRectangle(SVector(-2.0), SVector(2.0))
ABS.add_cells!(mng, UDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, UDom)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, labelset, UDom)

tstep = 1.0
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(1.0, u[1])
L_growthbound(u) = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)
sysnoise = SVector(1.0, 1.0)*0.001
measnoise = SVector(1.0, 1.0)*0.001

contsys = ABS.ControlSystemGrowthRK4(
    tstep,
    F_sys, L_growthbound,
    sysnoise, measnoise,
    nsys, ngrowthbound)

autom = ABS.AddAutomaton!(mng)
ABS.add_transitions!(mng, autom, stateset, labelset, stateset, contsys)

XInitDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9))
ABS.add_cells!(mng, XInitDom, rect, ABS.OUTER)
initset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, initset, XInitDom)

XTargetDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(0.0, 0.0), SVector(4.0, 4.0))
ABS.add_cells!(mng, XTargetDom, rect, ABS.OUTER)
targetset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, targetset, XTargetDom)

contr = ABS.AddController!(mng)
spec = ABS.ReachSpec(initset, targetset)
ABS.add_controls!(mng, contr, autom, spec)

xcell = ABS.get_some_cell(mng, XDom)
get_xcoord = ABS.cell2coord(mng.XDisc)
get_state = ABS.cell2symbol(mng.XMap)
x0 = get_xcoord(xcell)
source = get_state(xcell)
SLT_ = Any[]

for i in 1:6
    global source
    sourceset = ABS.AddStateSet!(mng)
    ABS.add_symbol!(mng, sourceset, source)
    labelset = ABS.AddLabelSet!(mng)
    ABS.add_symbols!(mng, labelset, contr, sourceset)
    targetset = ABS.AddStateSet!(mng)
    ABS.add_symbols!(mng, targetset, autom, sourceset, labelset)
    push!(SLT_, (sourceset, labelset, targetset))
    source = ABS.get_some_symbol(mng, targetset)
    source === nothing && break
end

Xs = ABS.AddXDomain!(mng)
Ys = ABS.AddXDomain!(mng)
Us = ABS.AddUDomain!(mng)

fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-5.5, 5.5))
ax.set_ylim((-5.3, 5.3))
Plot.cells!(ax, 1:2, mng, XDom, fa = 0.1)
Plot.cells!(ax, 1:2, mng, XInitDom)
Plot.cells!(ax, 1:2, mng, XTargetDom)
for (sourceset, labelset, targetset) in SLT_
    Plot.cells!(ax, 1:2, mng, sourceset, fc = "green")
    Plot.cells!(ax, 1:2, mng, targetset, fc = "blue")
    XDom = ABS.AddXDomain!(mng)
    ABS.add_cells!(mng, XDom, sourceset)
    UDom = ABS.AddUDomain!(mng)
    ABS.add_cells!(mng, UDom, labelset)
    Plot.images!(ax, 1:2, mng, XDom, UDom, contsys)
    # Same as above but to check with SymbolSets
    Plot.images!(ax, 1:2, mng, sourceset, labelset, contsys, fa = 0.0, ec = "red")
    # Plot.cell_approx!(ax, 1:2, Xs, Us, contsys)
end
Plot.trajectory!(ax, 1:2, mng, contr, contsys, x0, 10)

end  # module TestMain
