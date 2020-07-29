include("../src/abstraction.jl")

module TestMain

import Main.Abstraction
using PyPlot
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

lb = [0.0, 0.0]
ub = [10.0, 11.0]
x0 = [0.0, 0.0]
h = [1.0, 2.0]
X_grid = AB.NewGridSpaceHash(x0, h)
AB.add_to_gridspace!(X_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

lb = [-1.0]
ub = [1.0]
u0 = [0.0]
h = [0.5]
U_grid = AB.NewGridSpaceHash(u0, h)
AB.add_to_gridspace!(U_grid, AB.HyperRectangle(lb, ub), AB.OUTER)

trans_map = AB.NewTransitionMapHash(X_grid, U_grid, X_grid)
AB.add_transition_by_ref!(trans_map, 5, 6, 9)
AB.add_transition_by_ref!(trans_map, 5, 6, 10)
AB.add_transition_by_ref!(trans_map, 5, 7, 45)
AB.add_transition_by_ref!(trans_map, 15, 6, 45)
AB.add_transition_by_ref!(trans_map, 5, 6, 5)
AB.add_transition_by_ref!(trans_map, 15, 7, 45)
display(trans_map)

display(AB.get_transition_image(trans_map, 5))
display(AB.get_transition_image(trans_map, 15, 6))
display(AB.get_transition_image(trans_map, 15, 5))

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
