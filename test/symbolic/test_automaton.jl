module TestMain

using Test
using Dionysos
const DI = Dionysos
const SY = DI.Symbolic

sleep(0.1)
println("Started tests for automata")

function run_automaton_tests(AutomatonConstructor::Function)
    @testset "Automaton (impl = $(AutomatonConstructor))" begin
        nstates = 10
        nsymbols = 11
        autom = AutomatonConstructor(nstates, nsymbols)

        # Single transition per (state, input)
        SY.add_transition!(autom, 5, 9, 6)
        @test SY.is_deterministic(autom) == true

        # Now insert conflict on (5, 6)
        SY.add_transition!(autom, 5, 8, 6)
        @test SY.is_deterministic(autom) == false

        SY.add_transition!(autom, 5, 3, 7)
        SY.add_transition!(autom, 8, 3, 6)
        SY.add_transition!(autom, 5, 5, 6)
        SY.add_transition!(autom, 8, 3, 7)

        @test SY.ntransitions(autom) == 6

        SY.add_transitions!(autom, [(1, 2, 5), (1, 3, 4)])
        @test SY.ntransitions(autom) == 8

        targetlist = SY.post(autom, 5, 6)
        @test length(targetlist) == 3

        soursymblist = SY.pre(autom, 3)
        @test length(soursymblist) == 3
        @test collect(soursymblist) == [(5, 7), (8, 6), (8, 7)]

        @test length(SY.pre(autom, 4)) == 0
        @test collect(SY.pre(autom, 8)) == [(5, 6)]

        @test SY.is_deterministic(autom) == false
    end
end

# === Run tests for all known implementations ===

run_automaton_tests((n, m) -> SY.NewSortedAutomatonList(n, m))
run_automaton_tests((n, m) -> SY.NewIndexedAutomatonList(n, m))

println("End of tests")

end # module
