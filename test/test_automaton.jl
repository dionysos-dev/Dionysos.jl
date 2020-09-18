include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "Automaton" begin
nstates = 10
nsymbols = 11
autom = AB.NewAutomatonList(nstates, nsymbols)

AB.add_transition!(autom, 5, 6, 9)
AB.add_transition!(autom, 5, 6, 8)
AB.add_transition!(autom, 5, 7, 3)
AB.add_transition!(autom, 8, 6, 3)
AB.add_transition!(autom, 5, 6, 5)
AB.add_transition!(autom, 8, 7, 3)
@test AB.get_ntrans(autom) == 6
AB.add_transitions!(autom, [(1, 2, 5), (1, 3, 4)])
@test AB.get_ntrans(autom) == 8

targetlist = Int[]

AB.compute_post!(targetlist, autom, 5, 6)
@test length(targetlist) == 3
AB.compute_post!(targetlist, autom, 8, 6)
AB.compute_post!(targetlist, autom, 8, 5)
@test length(targetlist) == 4


soursymblist = AB.pre(autom, 3)
@test length(soursymblist) == 3
@test collect(soursymblist) == [(5, 7), (8, 6), (8, 7)]
soursymblist = AB.pre(autom, 4)
@test length(soursymblist) == 0
soursymblist = AB.pre(autom, 8)
@test length(soursymblist) == 1
@test collect(soursymblist)[1] == (5, 6)
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
