module TestMain

using Test
using Dionysos
const DI = Dionysos
const ST = DI.System

import HybridSystems

sleep(0.1)
println("Started tests for automata")

function run_automaton_tests(AutomatonConstructor::Function)
    @testset "Automaton (impl = $(AutomatonConstructor))" begin
        nstates = 10
        nsymbols = 11
        autom = AutomatonConstructor(nstates, nsymbols)

        # Single transition per (state, input)
        ST.add_transition!(autom, 5, 9, 6)
        @test ST.is_deterministic(autom) == true

        # Now insert conflict on (5, 6)
        ST.add_transition!(autom, 5, 8, 6)
        @test ST.is_deterministic(autom) == false

        ST.add_transition!(autom, 5, 3, 7)
        ST.add_transition!(autom, 8, 3, 6)
        ST.add_transition!(autom, 5, 5, 6)
        ST.add_transition!(autom, 8, 3, 7)

        @test HybridSystems.ntransitions(autom) == 6

        ST.add_transitions!(autom, [(1, 2, 5), (1, 3, 4)])
        @test HybridSystems.ntransitions(autom) == 8

        targetlist = ST.post(autom, 5, 6)
        @test length(targetlist) == 3

        soursymblist = ST.pre(autom, 3)
        @test length(soursymblist) == 3
        @test collect(soursymblist) == [(5, 7), (8, 6), (8, 7)]

        @test length(ST.pre(autom, 4)) == 0
        @test collect(ST.pre(autom, 8)) == [(5, 6)]

        @test ST.is_deterministic(autom) == false
    end
end

# === Run tests for all known implementations ===

run_automaton_tests((n, m) -> ST.NewSortedAutomatonList(n, m))
run_automaton_tests((n, m) -> ST.NewIndexedAutomatonList(n, m))

println("End of tests")

end # module
