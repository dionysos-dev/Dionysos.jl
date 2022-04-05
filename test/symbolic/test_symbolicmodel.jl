module TestMain

using Test
using StaticArrays
using ..Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "SymbolicModel" begin
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = DO.GridFree(x0, h)
Xfull = DO.DomainList(Xgrid)
DO.add_pos!(Xfull, (1, 1))
DO.add_pos!(Xfull, (2, 2))

u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = DO.GridFree(u0, h)
Ufull = DO.DomainList(Ugrid)
DO.add_pos!(Ufull, (0,))

symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
stateslist = Int[]
push!(stateslist, SY.get_state_by_xpos(symmodel, (1, 1)))
push!(stateslist, SY.get_state_by_xpos(symmodel, (2, 2)))
sort!(stateslist)
@test all(stateslist .== [1, 2])

xposlist = Tuple{Int, Int}[]
push!(xposlist, SY.get_xpos_by_state(symmodel, 1))
push!(xposlist, SY.get_xpos_by_state(symmodel, 2))
sort!(xposlist)
@test all(xposlist .== [(1, 1), (2, 2)])

uposlist = Tuple{Int}[]
push!(uposlist, SY.get_upos_by_symbol(symmodel, 1))
sort!(uposlist)
@test all(uposlist .== [(0,)])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
