include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "Controller" begin
contr = AB.NewControllerList()

AB.add_pair!(contr, 5, 6)
AB.add_pair!(contr, 5, 6)
AB.add_pair!(contr, 5, 7)
AB.add_pair!(contr, 8, 6)
AB.add_pair!(contr, 5, 6)
AB.add_pair!(contr, 8, 7)

symbollist = Int[]

AB.compute_enabled_symbols!(symbollist, contr, 5)
sort!(unique!(symbollist))
@test all(symbollist .== [6, 7])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
