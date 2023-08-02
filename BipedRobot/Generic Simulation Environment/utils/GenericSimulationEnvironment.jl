module GenericSimulationEnvironment
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

packagepath() = joinpath(@__DIR__, "../../ZMPBipedRobot/", "deps")

include("VirtualRobot.jl")
export VirtualRobot

end # End Module 
