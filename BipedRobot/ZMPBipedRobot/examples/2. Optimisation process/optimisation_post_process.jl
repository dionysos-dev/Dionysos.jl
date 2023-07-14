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

saveFolder = "examples/2. Optimisation process/results"

# read data from file into a matrix
# objectives = readdlm(saveFolder*"/objectives_values.txt", ',');
# solutions =  readdlm(saveFolder*"/solutions.txt", ',');
 
# Get the optimised paramters for a controller samplled at 50 Hz 
objectives = readdlm(saveFolder*"/objectives_values_Fs50.txt", ',')';
solutions =  readdlm(saveFolder*"/solutions_Fs50.txt", ',')';


# obj_values is a matrix of objective function values for each candidate
# where each row corresponds to a candidate and each column corresponds
# to an objective function
obj_values = objectives
n_candidates = size(obj_values, 1)

no_sorted_obj1 = Dict{Int, Vector{Float64}}() # create an empty dictionary

# Define comparison function based on first element of vector
compare(x, y) = x[1] < y[1]

for i = 1:n_candidates  # loop over the indices
    no_sorted_obj1[i] = obj_values[i, :] # store the objectives values for the index i in the dictionary
end

sorted_obj1 = sort(collect(no_sorted_obj1), by = x -> x[2][1])

obj1 = Float64[]
obj2 = Float64[]
obj3 = Float64[]
for (key, value) in sorted_obj1
    push!(obj1, value[1])
    push!(obj2, value[2])
    push!(obj3, value[3])
end
obj1 = obj1/maximum(obj1)
obj2 = obj2/maximum(obj2[isfinite.(obj2)])
obj3 = obj3/maximum(obj3[isfinite.(obj3)])

# obj1 = obj1[isfinite.(obj3)]/maximum(obj1)
# obj2 = obj2[isfinite.(obj3)]/maximum(obj2[isfinite.(obj2)])
# obj3 = obj3[isfinite.(obj3)]/maximum(obj3[isfinite.(obj3)])


plt = plot( #title = "Simualtion results for sorted by Walking Speed",
            xlabel = "Candidate", ylabel = "Normelized Objectives values",
            #xlim = best_range,
            legend=:top,  dpi=600)
msize = 1.75
scatter!(obj1, label = L"F_1(x)", ms = msize, msw = 0)
scatter!(obj2, label = L"F_2(x)", ms = msize, shape= :utriangle,  msw = 0)
scatter!(obj3, label = L"F_3(x)", ms = msize, shape = :rect, msw = 0)
display(plt)

# savefig(plt,saveFolder*"/Simualtion results for sorted by Walking Speed.png")

best_range = [1; 1512]
argmax(obj1[best_range[1] : best_range[2]])
best = argmin(obj2[best_range[1] : best_range[2]])
argmin(obj3[best_range[1] : best_range[2]])

best_objective = sorted_obj1[best_range[1] .+ best .- 1][2]
best_key = sorted_obj1[best_range[1] .+ best .- 1][1]

best_x = solutions[:, best_key]


plt = plot( #title = "Simualtion results for sorted by Walking Speed",
            xlabel = L"F_1(x)", ylabel = L"F_2(x)", zlabel = L"F_3(x)",
            #xlim = (1550,1570),
            #xlim = best_range,
            legend=:top,  dpi=600)
# scatter!(obj, label = L"F(x)", ms = msize, msw = 0)
scatter!(obj1,obj2,obj3, label = false, ms = msize, msw = .1)
display(plt) 

plt = plot( #title = "Simualtion results for sorted by Walking Speed",
            dpi=600,
            # layout = (3,1)
            )

# scatter!(plt[1], obj1, obj2, label = false, ms = msize, msw = .1, xlabel = L"F_1(x)", ylabel = L"F_2(x)",)
# scatter!(plt[2], obj2, obj3, label = false, ms = msize, msw = .1, xlabel = L"F_2(x)", ylabel = L"F_3(x)",)
scatter!(obj1, obj3, label = false, ms = msize, msw = .1, xlabel = L"F_1(x)", ylabel = L"F_3(x)",)
display(plt) 