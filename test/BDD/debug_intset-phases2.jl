include(joinpath(@__DIR__, "../../src/BDD/BDD.jl"))
using Test
using CUDD

println("")

set = BDD.IntSet()
BDD._phases2!(set, 3, 1000)
display(set.variables)
display(set.phases1)
display(set.phases2)
display(set.newphases1)
display(set.newphases2)
display(set.auxphases)
display(set.newvars)

print("")
