module Utils

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Polyhedra
import LazySets
import SpecialFunctions
using LinearAlgebra

include("plotting/simple_plots.jl")

include("functions.jl")
include("periodic.jl")
include("incl_mode.jl")
include("data_structures/sorted_vector_set.jl")
include("data_structures/tree.jl")
include("search/RRT.jl")

include("optim/scalar_optimization.jl")

include("sets/lazy_set_operations.jl")
include("sets/rectangle.jl")
include("sets/ellipsoid.jl")
include("sets/degenerate_ellipsoid.jl")
include("sets/semilinear_set.jl")

include("pclf.jl")

end  # module Utils
