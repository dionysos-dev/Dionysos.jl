module Utils

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Colors
import Polyhedra
import DataStructures
import LazySets
import IntervalArithmetic as IA
import SpecialFunctions: gamma
using LinearAlgebra, JuMP
import PlotUtils

include("plotting/colorbar.jl")
include("plotting/simple_plots.jl")

include("scalar_functions.jl")
include("data_structures/sorted_vector_set.jl")
include("data_structures/queue.jl")
include("data_structures/tree.jl")
include("data_structures/digraph.jl")
include("search/generic_search.jl")
include("search/RRT.jl")

include("optim/branch_and_bound.jl")
include("optim/bisection.jl")
include("optim/newton_method.jl")

include("sets/lazy_set_operations.jl")
include("sets/rectangle.jl")
include("sets/ellipsoid.jl")
include("sets/degenerate_ellipsoid.jl")
include("sets/polyhedron.jl")

include("box.jl")
using JuMP
include("ellipsoidal_transitions.jl")

include("pclf.jl")

end  # module Utils
