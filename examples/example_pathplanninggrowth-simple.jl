include("pathplanning.jl")
using ProfileView
@profview PathPlanning.path_planning(4.0, nstep = 100, approx_mode = "growth")
