## Import useful packages 
using LinearAlgebra
using StaticArrays
using StructArrays
using Plots
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using StaticArrays
using Symbolics
using MeshCat, MeshCatMechanisms, Blink
using MechanismGeometries
using LaTeXStrings
using Tables, CSV
using LightXML
using DelimitedFiles, DataFrames

## Include and import the ZMP based controller 
include(joinpath(@__DIR__, "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot

packagepath() = ZMProbot.packagepath() 
urdfpath() = joinpath(packagepath(), "ZMP_2DBipedRobot.urdf")
# include(joinpath(packagepath(), "param.jl"))

visual = URDFVisuals(urdfpath())
e = root(visual.xdoc)

temp = e["joint"][5]
temp["parent"]

println(attribute(temp, "name"))
        
L1 = 0.0 
L2 = 0.0 
d = 0.0 
offset_hip_to_motor = 0.0 
offset_ankle_to_foot = 0.0

for  (i, joint) in enumerate(e["joint"]) 
    if (attribute(joint, "name")) == "l_hip_to_motor"
        println(i)
        s = attribute(joint["origin"]..., "xyz")
        numbers = split(s, " ")
        offset_hip_to_motor = abs(parse(Float64, numbers[3])) 
        d = 2 * abs(parse(Float64, numbers[2])) 
    end 
    if  (attribute(joint, "name")) == "l_thigh_link_to_motor"
        println(i)
        s = attribute(joint["origin"]..., "xyz")
        numbers = split(s, " ")
        L1 = abs(parse(Float64, numbers[3])) 
    end 
    if  (attribute(joint, "name")) == "l_ankle"
        println(i)
        s = attribute(joint["origin"]..., "xyz")
        numbers = split(s, " ")
        L2 = abs(parse(Float64, numbers[3])) 
    end 
end 

for (i, link) in enumerate(e["link"])
    if (attribute(link, "name") == "l_foot_link")
        println(i)
        s = attribute(link["visual"][1]["geometry"][1]["box"]..., "size")
        numbers = split(s, " ")
        offset_ankle_to_foot = abs(parse(Float64, numbers[3])) 
    end 
    if (attribute(link, "name") == "base_boom_link")
        println(i)
        s = attribute(link["visual"][1]["geometry"][1]["cylinder"]..., "size")
        # numbers = split(s, " ")
    end 
end 

println(d)
println(L1)
println(L2)
println(offset_hip_to_motor)
println(offset_ankle_to_foot)