module TestMain

using Test
using Dionysos
const DI = Dionysos
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "Proba Automaton" begin
    nstates = 10
    nsymbols = 11
    autom = SY.NewProbaAutomaton(nstates, nsymbols)

    SY.add_transition!(autom, 5, 9, 6, 0.5)
    SY.add_transition!(autom, 5, 8, 6, 0.1)
    SY.add_transition!(autom, 5, 3, 7, 1.0)
    SY.add_transition!(autom, 8, 3, 6, 0.2)
    SY.add_transition!(autom, 5, 5, 6, 0.3)
    SY.add_transition!(autom, 8, 3, 7, 0.8)
    @test SY.ntransitions(autom) == 6
    SY.add_transitions!(autom, [(1, 2, 5, 0.5), (1, 3, 4, 0.1)])
    @test SY.ntransitions(autom) == 8

    soursymblist = SY.pre(autom, 6)
    @test length(soursymblist) == 4
    @test collect(soursymblist) == [(5, 5, 0.3), (5, 8, 0.1), (5, 9, 0.5), (8, 3, 0.2)]
    soursymblist = SY.pre(autom, 4)
    @test length(soursymblist) == 0
    soursymblist = SY.pre(autom, 8)
    @test length(soursymblist) == 0

    SY.empty!(autom)
    @test SY.ntransitions(autom) == 0
end

println("End test")

end  # module TestMain
