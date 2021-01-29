include("../../src/Abstraction/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "Controller" begin
contr = AB.NewControllerList()

AB.push_new!(contr, (5, 6))
AB.push_new!(contr, (5, 6))
AB.push_new!(contr, (5, 7))
AB.push_new!(contr, (8, 6))
AB.push_new!(contr, (5, 6))
AB.push_new!(contr, (8, 7))

symbollist = [x[1] for x in AB.fix_and_eliminate_first(contr, 5)]
sort!(unique!(symbollist))
@test all(symbollist .== [6, 7])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
