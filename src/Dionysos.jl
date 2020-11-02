module Dionysos

# TODO Reenable once https://github.com/sisl/CUDD.jl/issues/15#issuecomment-719958808 is resolved
#include("BDD/BDD.jl")

include("optimal_control.jl")
include("bemporad_morari.jl")
include("branch_and_bound.jl")

end # module
