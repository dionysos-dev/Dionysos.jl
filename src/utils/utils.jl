module Utils

using StaticArrays
using LinearAlgebra


include("files/files_management.jl")
##### PLOTS
using Plots, Colors
include("plotting/colorbar.jl")
include("plotting/simple_plots.jl")
#####

# Temporarily remove BDDs
# include("data_structures/BDD/BDD.jl")
include("data_structures/sorted_vector_set.jl")
include("data_structures/queue.jl")
include("data_structures/tree.jl")
include("data_structures/digraph.jl")
include("search/generic_search.jl")
include("search/RRT.jl")


include("optim/branch_and_bound.jl")
include("optim/bisection.jl")
include("optim/newton_method.jl")

include("rectangle.jl")
include("box.jl")
include("ellipsoid.jl")
include("degenerate_ellipsoid.jl")
include("polyhedron.jl")
include("monteCarlo.jl")

include("lazySetOperations.jl")

end  # module Utils
