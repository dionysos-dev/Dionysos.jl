include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "SymbolicModel" begin
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = AB.GridFree(x0, h)
Xfull = AB.DomainList(Xgrid)
AB.add_pos!(Xfull, (1, 1))
AB.add_pos!(Xfull, (2, 2))

u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = AB.GridFree(u0, h)
Ufull = AB.DomainList(Ugrid)
AB.add_pos!(Ufull, (0,))

symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)
stateslist = Int[]
push!(stateslist, AB.get_state_by_xpos(symmodel, (1, 1)))
push!(stateslist, AB.get_state_by_xpos(symmodel, (2, 2)))
sort!(stateslist)
@test all(stateslist .== [1, 2])

xposlist = Tuple{Int, Int}[]
push!(xposlist, AB.get_xpos_by_state(symmodel, 1))
push!(xposlist, AB.get_xpos_by_state(symmodel, 2))
sort!(xposlist)
@test all(xposlist .== [(1, 1), (2, 2)])

uposlist = Tuple{Int}[]
push!(uposlist, AB.get_upos_by_symbol(symmodel, 1))
sort!(uposlist)
@test all(uposlist .== [(0,)])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
