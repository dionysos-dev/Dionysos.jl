module Utils

using StaticArrays

include("data_structures/BDD/BDD.jl")
include("data_structures/sorted_vector_set.jl")

include("data_structures/queue.jl")
include("search/generic_search.jl")

include("optim/branch_and_bound.jl")

include("rectangle.jl")
include("polyhedron.jl")


end  # module Utils
