using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
using LinearAlgebra
using JLD2
import MathOptInterface as MOI

include("../robot_vcontrol.jl")
import .RobotVelocity as RV

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

function save_optimizer(filename::AbstractString, optimizer)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        return file["optimizer"] = optimizer
    end
    return filename
end

function load_optimizer(filename::AbstractString)
    return JLD2.jldopen(filename, "r") do file
        return file["optimizer"]
    end
end

# load_optimizer(joinpath(@__DIR__, "out/optimizer.jld2"))

filename = joinpath(@__DIR__, "out/optimizer.jld2")
optimizer = JLD2.load(filename, "optimizer")

MOI.optimize!(optimizer)

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
println("OptimalControlProblem success: $success")

if success
    # save_optimizer(joinpath(@__DIR__, "out/optimizer.jld2"), optimizer)
    JLD2.jldsave(filename; optimizer = optimizer)
end