import Spot

# 1) Define an LTL formula (atomic props are symbols like :goal, :unsafe, etc.)
φ = Spot.ltl"G(!unsafe) & F(goal)"
println("Formula: ", φ)

# 2) Build a deterministic Rabin automaton (pure Julia structure in Spot.jl)
dra = Spot.DeterministicRabinAutomata(φ)
println("Built a DRA: ", typeof(dra))

# 3) Test stepping the automaton with different valuations
#    Valuation format: a tuple of symbols that are TRUE at the current step.
#    Example: (:goal,) means goal is true, unsafe is false.
q0 = 1  # initial automaton state index (Spot.jl examples use small integers)
println("Initial DRA state q0 = ", q0)

vals = [
    (),              # nothing true
    (:unsafe,),      # unsafe true
    (:goal,),        # goal true
    (:goal, :unsafe) # both true
]

println("\nNext-state tests:")
for v in vals
    q1 = Spot.nextstate(dra, q0, v)
    println("  δ($q0, $v) = $q1")
end

