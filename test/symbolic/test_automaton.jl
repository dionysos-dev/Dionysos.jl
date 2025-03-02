module TestMain

using Test
using Dionysos
const DI = Dionysos
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

sleep(0.1) # used for good printing
println("Started test")

@testset "Automaton" begin
    nstates = 10
    nsymbols = 11
    autom = SY.NewAutomatonList(nstates, nsymbols)

    SY.add_transition!(autom, 5, 9, 6)
    @test SY.is_deterministic(autom) == true
    SY.add_transition!(autom, 5, 8, 6)
    @test SY.is_deterministic(autom) == false
    SY.add_transition!(autom, 5, 3, 7)
    SY.add_transition!(autom, 8, 3, 6)
    SY.add_transition!(autom, 5, 5, 6)
    SY.add_transition!(autom, 8, 3, 7)
    @test SY.ntransitions(autom) == 6
    SY.add_transitions!(autom, [(1, 2, 5), (1, 3, 4)])
    @test SY.ntransitions(autom) == 8

    targetlist = Int[]

    SY.compute_post!(targetlist, autom, 5, 6)
    @test length(targetlist) == 3
    SY.compute_post!(targetlist, autom, 8, 6)
    SY.compute_post!(targetlist, autom, 8, 5)
    @test length(targetlist) == 4

    soursymblist = SY.pre(autom, 3)
    @test length(soursymblist) == 3
    @test collect(soursymblist) == [(5, 7), (8, 6), (8, 7)]
    soursymblist = SY.pre(autom, 4)
    @test length(soursymblist) == 0
    soursymblist = SY.pre(autom, 8)
    @test length(soursymblist) == 1
    @test collect(soursymblist)[1] == (5, 6)

    @test SY.is_deterministic(autom) == false
end

println("End test")

end  # module TestMain
