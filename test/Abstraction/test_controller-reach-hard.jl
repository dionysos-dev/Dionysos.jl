include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/controller-reach" begin
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
@test ABS.get_ntransitions(mng, autom) === 15530

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
spec = ABS.ReachSpec(mng, initset, targetset)
@test ABS.add_controls!(mng, contr, autom, spec)
@test ABS.get_npairs(mng, contr) === 412

if VERSION >= v"1.5"
    function ftestallocate(mng, contr, autom, initset, targetset)
        npostsUNDRV, ninitsDRV, allDRV, currentDRV, nextDRV =
            ABS._initialize_reach(mng, autom, initset, targetset)
        # Preallocates to make sure `_compute_controller_reach` does not need to allocate
        sizehint!(contr.pairs, 10_000)
        sizehint!(allDRV, 10_000)
        sizehint!(currentDRV, 10_000)
        sizehint!(nextDRV, 10_000)
        return @allocated ABS._add_controls_reach!(mng, contr, autom, initset,
            npostsUNDRV, ninitsDRV, allDRV, currentDRV, nextDRV)
    end
    autom = ABS.AddAutomaton!(mng)
    ABS.add_transitions!(mng, autom, stateset, labelset, stateset)
    contr = ABS.AddController!(mng)
    ftestallocate(mng, contr, autom, stateset, targetset)
    contr = ABS.AddController!(mng)
    @test ftestallocate(mng, contr, autom, stateset, targetset) === 0
end
print("")
end

end  # module TestMain
