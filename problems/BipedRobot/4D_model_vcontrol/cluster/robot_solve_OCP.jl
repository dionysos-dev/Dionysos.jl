import Pkg
Pkg.instantiate()
 
using Printf

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

# function save_optimizer(filename::AbstractString, optimizer)
#     mkpath(dirname(filename))
#     JLD2.jldopen(filename, "w") do file
#         return file["optimizer"] = optimizer
#     end
#     return filename
# end

# function load_optimizer(filename::AbstractString)
#     return JLD2.jldopen(filename, "r") do file
#         return file["optimizer"]
#     end
# end

# -----------------------
# Load optimizer
# -----------------------

outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", "")
optimizer = JLD2.load(outfile, "optimizer")

# ------------------------------------------------------------------------------
# Problem (OCP)
# ------------------------------------------------------------------------------

include(joinpath(@__DIR__, "robot_create_problem.jl"))

concrete_problem =
    DI.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        0.0
    )

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

# ------------------------------------------------------------------------------
# Solve OCP
# ------------------------------------------------------------------------------

@info "Solving the OCP on the loaded abstraction"

t_solve = @elapsed MOI.optimize!(optimizer)

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

if success
    t_save = @elapsed JLD2.jldsave(outfile; optimizer = optimizer)
    @printf("OptimalControlProblem success!")
    @printf("Solve time:                   %.3f s\n", t_solve)
    @printf("Save time:                    %.3f s\n", t_save)
    @info "Saved solved optimizer" outfile
else
    @printf("No solution found for OptimalControlProblem")
end