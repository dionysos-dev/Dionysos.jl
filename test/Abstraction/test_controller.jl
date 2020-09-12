include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/controller" begin
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

contr = ABS.AddController!(mng)
CoT = ABS.controltype(mng)
@test ABS.get_npairs(mng, contr) === 0
ABS.add_controls!(mng, contr, stateset, labelset)
nstates = ABS.get_nstates(mng)
nlabels = ABS.get_nlabels(mng)
@test ABS.get_npairs(mng, contr) === nstates*nlabels
pairlist = Set(CoT(state, label)
    for state in ABS.enum_symbols(mng, stateset),
        label in ABS.enum_symbols(mng, labelset))
@test Set(ABS.enum_pairs(mng, contr)) == pairlist
enablelabelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, enablelabelset, contr, stateset)
@test Set(ABS.enum_symbols(mng, labelset)) == Set(ABS.enum_symbols(mng, enablelabelset))
state = ABS.get_some_symbol(mng, stateset)
@test Set(ABS.enum_symbols(mng, labelset)) == Set(ABS.enum_enabled_labels(mng, contr, state))
# pair = first(contr.pairs)
# Main.@code_warntype hash(pair, UInt(0))
@test empty!(mng, contr) === contr
@test ABS.get_npairs(mng, contr) === 0
print("")
end

end  # module TestMain
