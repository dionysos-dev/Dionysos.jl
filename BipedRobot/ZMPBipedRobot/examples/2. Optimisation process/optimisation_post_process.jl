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

using DelimitedFiles

local_dir = joinpath(@__DIR__, "..", "../")
saveFolder = local_dir * "docs/2. Optimisation process/"

# Get the optimised paramters for a controller samplled at 50 Hz 
objectives = readdlm(saveFolder * "/objectives_values_Fs50.txt", ',')';
solutions = readdlm(saveFolder * "/solutions_Fs50.txt", ',')';

# obj_values is a matrix of objective function values for each candidate
# where each row corresponds to a candidate and each column corresponds
# to an objective function
obj_values = objectives
n_candidates = size(obj_values, 1)

no_sorted_obj1 = Dict{Int, Vector{Float64}}() # create an empty dictionary

for i in 1:n_candidates  # loop over the indices
    no_sorted_obj1[i] = obj_values[i, :] # store the objectives values for the index i in the dictionary
end

sorted_obj1 = sort(collect(no_sorted_obj1); by = x -> x[2][1])

obj1 = Float64[]
obj2 = Float64[]
obj3 = Float64[]
for (key, value) in sorted_obj1
    push!(obj1, value[1])
    push!(obj2, value[2])
    push!(obj3, value[3])
end
obj1 = obj1 / maximum(obj1)
obj2 = obj2 / maximum(obj2[isfinite.(obj2)])
obj3 = obj3 / maximum(obj3[isfinite.(obj3)])
best_range = [1, length(obj1)]
plt = plot(;
    xlabel = "Candidate",
    ylabel = "Normelized Objectives values",
    xlim = best_range,
    legend = :top,
    dpi = 600,
)
msize = 1.75
scatter!(obj1; label = L"F_1(x)", ms = msize, msw = 0)
scatter!(obj2; label = L"F_2(x)", ms = msize, shape = :utriangle, msw = 0)
scatter!(obj3; label = L"F_3(x)", ms = msize, shape = :rect, msw = 0)
display(plt)

best_range = [1; 800]
argmax(obj1[best_range[1]:best_range[2]])
best = argmin(obj2[best_range[1]:best_range[2]])
argmin(obj3[best_range[1]:best_range[2]])

best_objective = sorted_obj1[best_range[1] .+ best .- 1][2]
best_key = sorted_obj1[best_range[1] .+ best .- 1][1]

best_x = solutions[:, best_key]
