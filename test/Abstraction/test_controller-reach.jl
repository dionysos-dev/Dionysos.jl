include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/controller-reach" begin
xorig = SVector(0)
xh = SVector(1)
uorig = SVector(0)
uh = SVector(1)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

XDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(0), SVector(3))
ABS.add_cells!(mng, XDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, XDom)
stateset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, stateset, XDom)

UDom = ABS.AddUDomain!(mng)
rect = ABS.HyperRectangle(SVector(0), SVector(3))
ABS.add_cells!(mng, UDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, UDom)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, labelset, UDom)

autom = ABS.AddAutomaton!(mng)
ABS.add_transitions!(mng, autom, stateset, labelset, stateset)

contr = ABS.AddController!(mng)
spec = ABS.ReachSpec(mng, stateset, stateset)
@test ABS.add_controls!(mng, contr, autom, spec)
@test ABS.get_npairs(mng, contr) === 0

autom = ABS.AddAutomaton!(mng)
targetset = ABS.AddStateSet!(mng)
ABS.add_symbol!(mng, targetset, ABS.get_some_symbol(mng, stateset))
ABS.add_transitions!(mng, autom, stateset, labelset, targetset)
contr = ABS.AddController!(mng)
spec = ABS.ReachSpec(mng, stateset, targetset)
@test ABS.add_controls!(mng, contr, autom, spec)
@test ABS.get_npairs(mng, contr) === 3

autom = ABS.AddAutomaton!(mng)
statelist = collect(ABS.enum_symbols(mng, stateset))
labellist = collect(ABS.enum_symbols(mng, labelset))
TT = ABS.transitiontype(mng)
for i = 1:length(statelist)-1
    j = mod(i + 20, length(labellist)) + 1
    ABS.add_transition!(mng, autom, TT(statelist[i], labellist[j], statelist[i+1]))
end
targetset = ABS.AddStateSet!(mng)
ABS.add_symbol!(mng, targetset, statelist[end])
for i = 1:length(statelist)-1
    initset = ABS.AddStateSet!(mng)
    ABS.add_symbol!(mng, initset, statelist[i])
    contr = ABS.AddController!(mng)
    spec = ABS.ReachSpec(mng, initset, targetset)
    @test ABS.add_controls!(mng, contr, autom, spec)
    @test ABS.get_npairs(mng, contr) === length(statelist) - i
end

if VERSION >= v"1.5"
    function ftestallocate(mng, contr, autom, initset, targetset)
        npostsUNDRV, ninitsDRV, allDRV, currentDRV, nextDRV =
            ABS._initialize_reach(mng, autom, initset, targetset)
        # Preallocates to make sure `_compute_controller_reach` does not need to allocate
        sizehint!(contr.pairs, 1_000)
        sizehint!(allDRV, 1_000)
        sizehint!(currentDRV, 1_000)
        sizehint!(nextDRV, 1_000)
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
