include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/automaton" begin
xorig = SVector(0.0, 0.0)
xh = SVector(1.0, 2.0)
uorig = SVector(0.0)
uh = SVector(0.5)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

XDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(0.0, 0.0), SVector(10.0, 11.0))
ABS.add_cells!(mng, XDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, XDom)
stateset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, stateset, XDom)

UDom = ABS.AddUDomain!(mng)
rect = rect = ABS.HyperRectangle(SVector(-1.0), SVector(1.0))
ABS.add_cells!(mng, UDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, UDom)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, labelset, UDom)

autom = ABS.AddAutomaton!(mng)
TT = ABS.transitiontype(mng)
@test ABS.get_ntransitions(mng, autom) === 0
state = ABS.get_some_symbol(mng, stateset)
label = ABS.get_some_symbol(mng, labelset)
ABS.add_transition!(mng, autom, TT(state, label, state))
@test ABS.get_ntransitions(mng, autom) === 1
@test empty!(mng, autom) === autom
@test ABS.get_ntransitions(mng, autom) === 0
ABS.add_transitions!(mng, autom, stateset, labelset, stateset)
nstates = ABS.get_nstates(mng)
nlabels = ABS.get_nlabels(mng)
@test ABS.get_ntransitions(mng, autom) === nstates*nlabels*nstates
translist = Set(TT(source, label, target)
    for source in ABS.enum_symbols(mng, stateset),
        label in ABS.enum_symbols(mng, labelset),
        target in ABS.enum_symbols(mng, stateset))
@test Set(trans for trans in ABS.enum_transitions(mng, autom)) == translist
targetset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, targetset, autom, stateset, labelset)
@test Set(ABS.enum_symbols(mng, stateset)) == Set(ABS.enum_symbols(mng, targetset))
print("")
end

end  # module TestMain
