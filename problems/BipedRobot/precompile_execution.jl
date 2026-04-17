using Dionysos
using StaticArrays

include(joinpath(@__DIR__, "6D_model", "robot_problem.jl"))
using .RobotProblem

robot_urdf = joinpath(@__DIR__, "deps", "ZMP_2DBipedRobot_nodamping.urdf")

RobotProblem.warmup_robot_problem!(; robot_urdf = robot_urdf, tstep = 0.1)

sys = RobotProblem.system(; robot_urdf = robot_urdf, tstep = 0.1)
x = SVector{6, Float64}(zeros(6))
u = SVector{3, Float64}(zeros(3))
sys.f(x, u)
