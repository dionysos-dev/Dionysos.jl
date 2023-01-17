module Utils

using StaticArrays
using LinearAlgebra

include("data_structures/BDD/BDD.jl")
include("data_structures/sorted_vector_set.jl")

include("data_structures/queue.jl")
include("search/generic_search.jl")

include("optim/branch_and_bound.jl")
include("optim/bisection.jl")
include("optim/newton_method.jl")

include("rectangle.jl")
include("ellipsoid.jl")
include("polyhedron.jl")
include("monteCarlo.jl")

include("lazySetOperations.jl")

end  # module Utils
