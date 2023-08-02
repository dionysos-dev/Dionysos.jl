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

## Include and import the ZMP based controller 
include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
include(joinpath(@__DIR__, "..", "..", "deps", "param.jl"))
L1 = 0.20125
L2 = 0.172

import .ZMPBipedRobot
const ZMProbot = ZMPBipedRobot

local_dir = joinpath(@__DIR__, "..", "../")
saveFolder = local_dir * "docs/2. Optimisation process"

Δz_space = (0.05:0.05:0.3) * zmax
Tstep_space = 0.1:0.1:1
Lmax_space = 0.05:0.05:((L1 + L2) / 2)
δ_space = 0.05:0.05:0.2
δver_space = 0
hstep_space = 0.02:0.01:0.05
i = 1;
solutions = Array[]
Niterations =
    length(Δz_space) *
    length(Tstep_space) *
    length(Lmax_space) *
    length(δ_space) *
    length(δver_space) *
    length(hstep_space)

objectives = Array[]

for Δz in Δz_space
    for Tstep in Tstep_space
        for Lmax in Lmax_space
            for δ in δ_space
                for δver in δver_space
                    Tver = δver * Tstep
                    for hstep in hstep_space
                        x = [Δz Tstep Lmax δ Tver hstep]
                        if (x[3] <= sqrt(-(x[1]) * (x[1] - L1 - L2))) &&
                           (x[2] * x[4] > 2 * Ts)
                            obj, _, _ = ZMProbot.ZMPbipedObjective(x)
                            push!(objectives, obj)
                            push!(solutions, x)
                            println(
                                "+====================================================+",
                            )
                            println("Iteration : $(i)/$(Niterations)")
                            println("x = $(x)")
                            println("Objectives : $(obj)")
                            global i = i + 1
                        end
                    end
                end
            end
        end
    end
end
objectives = reduce(hcat, objectives)
plt_objective = plot(;
    xlabel = "max $(L"f_1(x)")",
    ylabel = "min $(L"f_2(x)")",
    zlabel = "min $(L"f_3(x)")",
    camera = (45, 45),
    dpi = 600,
)
scatter!(
    objectives[1, :],
    objectives[2, :],
    [objectives[3, :]];
    label = false,
    ms = 2,
    msw = 1,
)
display(plt_objective)

## Write into text file 
using DelimitedFiles
writedlm(saveFolder * "/objectives_values_Fs50.txt", objectives, ',')
writedlm(saveFolder * "/solutions_Fs50.txt", solutions, ',')
