module Utils

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Polyhedra
import LazySets
import SpecialFunctions
using LinearAlgebra

include("plotting.jl")

include("functions.jl")
include("periodic.jl")
include("incl_mode.jl")
include("data_structures/sorted_vector_set.jl")
include("data_structures/tree.jl")
include("search/RRT.jl")

include("numeric/scalar_optimization.jl")

include("sets/set_algebra.jl")
include("sets/ellipsoid.jl")
include("sets/semilinear_set.jl")

include("pclf.jl")

end  # module Utils
