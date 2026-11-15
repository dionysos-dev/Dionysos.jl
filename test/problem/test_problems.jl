module TestMain

using Test
using Dionysos
const DI = Dionysos
const PR = DI.Problem

sleep(0.1) # used for good printing
println("Started test")

@testset "Problems" begin
    @test Base.isfinite(4) == true
    @test Base.isfinite(-2.5) == true
    @test Base.isfinite(Inf) == false
end

println("End test")

end  # module TestMain
