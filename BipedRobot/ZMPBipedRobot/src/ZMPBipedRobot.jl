module ZMPBipedRobot
using LinearAlgebra
using StructArrays
using Plots
using CSV, Tables
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using StaticArrays
using Symbolics
using MeshCat, MeshCatMechanisms, Blink
using MechanismGeometries
using LightXML
using GeometryTypes
using Random
using DataStructures
using LaTeXStrings
using DataFrames

packagepath() = joinpath(@__DIR__, "..", "deps")

#     function __init__()
#         if !isfile(urdfpath())
#             error("Could not find $(urdfpath()). Please run `import Pkg; Pkg.build(\"ZMPBipedRobot\")`.")
#         end
#     end

include("BipedRobot.jl")
export BipedRobot

include("util.jl")
export cartTableModel, continuous2discrete, eye, getSplineCoeff, spline

include("PreviewControl.jl")
export PreviewController, computeGains

include("FootsPlacement.jl")
export FootPlanner, computeFootsPlacement

include("ZMPTrajectoryGenerator.jl")
export ZMPTrajectory, computeZMPTrajectory

include("CoMTrajectoryGenerator.jl")
export CoMTrajectory, computeCoMTrajectory

include("SwingFootTrajectoryGenerator.jl")
export SwingFootTrajectory, computeSwingFootTrajectory

include("DirectKinematics.jl")
export foot2hip

include("InverseKinematics.jl")
export InverseKinematics, computeGlobal2Local

include("RobotSimulator.jl")
export RobotSimulator,
    getMechanism,
    set_nominal!,
    set_initialbody!,
    update_visulizer!,
    show_frame!,
    trajectory_controller!,
    simulate,
    measureZMP
include("WalkingOptimization.jl")
export computeEnergy, simulate

include("positionControl.jl")
export pid_control!

end # module ZMPBipedRobot
