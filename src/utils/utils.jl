module Utils

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Polyhedra
import LazySets
import IntervalArithmetic as IA
using LinearAlgebra, JuMP

include("plotting/simple_plots.jl")

include("scalar_functions.jl")
include("data_structures/sorted_vector_set.jl")
include("data_structures/tree.jl")
include("search/RRT.jl")

include("optim/bisection.jl")
include("optim/newton_method.jl")

include("sets/lazy_set_operations.jl")
include("sets/rectangle.jl")
include("sets/ellipsoid.jl")
include("sets/degenerate_ellipsoid.jl")
include("sets/polyhedron.jl")

include("box.jl")
include("ellipsoidal_transitions.jl")

include("pclf.jl")

end  # module Utils
