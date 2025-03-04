module Utils

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Colors
import Polyhedra
import DataStructures
import LazySets
import IntervalArithmetic: IntervalBox
import SpecialFunctions: gamma
using LinearAlgebra, JuMP
import Plots: palette, annotate!
import HiGHS

include("files/files_management.jl")

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

include("rectangle.jl")
include("box.jl")
include("ellipsoid.jl")
include("degenerate_ellipsoid.jl")
include("polyhedron.jl")
include("intersection_set.jl")

include("lazy_set_operations.jl")

end  # module Utils
