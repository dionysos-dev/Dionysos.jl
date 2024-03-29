module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const OP = DI.Optim
const AB = OP.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "Controller" begin
    contr = AB.UniformGridAbstraction.NewControllerList()

    UT.push_new!(contr, (5, 6))
    UT.push_new!(contr, (5, 6))
    UT.push_new!(contr, (5, 7))
    UT.push_new!(contr, (8, 6))
    UT.push_new!(contr, (5, 6))
    UT.push_new!(contr, (8, 7))

    symbollist = [x[1] for x in UT.fix_and_eliminate_first(contr, 5)]
    sort!(unique!(symbollist))
    @test all(symbollist .== [6, 7])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
