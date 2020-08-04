include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "SymbolicModel" begin
x0 = (0.0, 0.0)
h = (1.0, 2.0)
Xgrid = AB.NewGridSpaceList(x0, h)
AB.add_pos!(Xgrid, (1, 1))
AB.add_pos!(Xgrid, (2, 2))

u0 = (0.0,)
h = (0.5,)
Ugrid = AB.NewGridSpaceList(u0, h)
AB.add_pos!(Ugrid, (0,))

symmodel = AB.NewSymbolicModelListList(Xgrid, Ugrid)
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
