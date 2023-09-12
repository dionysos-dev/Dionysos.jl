# # Ellipsoids
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Ellipsoids.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Ellipsoids.ipynb)
#
# In this file, we will visit the basic features provided by Dionysos for manipulating ellipsoids:
# - inclusion
# - intersection
# 
# First, let us import a few packages that are necessary to run this example.
using Dionysos
using StaticArrays
using LinearAlgebra
using Plots

# The main package [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) provides most important data structures that we will need.

const DI = Dionysos
const UT = DI.Utils

# We define a plotting functions
function plot_config!(fig, El0, El, Elnew)
    if Elnew ∈ El0
        plot!(fig, El0; color = :red, label = "El0")
        plot!(fig, Elnew; color = :green, label = "Elnew")
    else
        plot!(fig, Elnew; color = :green, label = "Elnew")
        plot!(fig, El0; color = :red, label = "El0")
    end
    return plot!(fig, El; color = :blue, label = "El", show = true)
end

function analyze(fig1, fig2, i)
    El0 = E0L[i]
    Elnew = UT.scale_for_inclusion_contact_point(El0, El)
    plot_config!(fig1, El0, El, Elnew)
    println(El0 ∈ El ? "El0 ∈ El" : "El0 ∉ El")

    Elnew = UT.scale_for_noninclusion_contact_point(El0, El)
    plot_config!(fig2, El0, El, Elnew)
    return println(UT.is_intersected(El0, El) ? "El0 ∩ El ≠ ∅" : "El0 ∩ El = ∅")
end

# We define some ellipsoids
c = [1.5; 1.5]
P = [
    4.0 0.5
    0.5 6.0
]
El = UT.Ellipsoid(P, c)

P0 = [
    0.4 -0.1
    -0.1 0.5
]
vals = [4.1, 3.32, 2.8, 2.4]
E0L = [UT.Ellipsoid(P0, [c0x; c0x - 0.2]) for c0x in vals]

# ### Case 1: non intersection
fig1_1 = plot(; aspect_ratio = :equal);
fig1_2 = plot(; aspect_ratio = :equal);
analyze(fig1_1, fig1_2, 1)

#
plot!(fig1_1)

#
plot!(fig1_2)

# ### Case 2: non intersection
fig2_1 = plot(; aspect_ratio = :equal);
fig2_2 = plot(; aspect_ratio = :equal);
analyze(fig2_1, fig2_2, 2)

#
plot!(fig2_1)

#
plot!(fig2_2)

# ### Case 3: intersection, non inclusion
fig3_1 = plot(; aspect_ratio = :equal);
fig3_2 = plot(; aspect_ratio = :equal);
analyze(fig3_1, fig3_2, 3)

#
plot!(fig3_1)

#
plot!(fig3_2)

# ### Case 4: inclusion
fig4_1 = plot(; aspect_ratio = :equal);
fig4_2 = plot(; aspect_ratio = :equal);
analyze(fig4_1, fig4_2, 4)

#
plot!(fig4_1)

#
plot!(fig4_2)
