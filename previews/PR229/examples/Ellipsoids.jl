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
function plot_config!(El0, El, Elnew)
    fig = plot(aspect_ratio=:equal);
    if Elnew∈El0
        plot!(fig, El0, color=:red, label="El0");
        plot!(fig, Elnew, color=:green, label="Elnew");
    else
        plot!(fig, Elnew, color=:green, label="Elnew");
        plot!(fig, El0, color=:red, label="El0");
    end
    plot!(fig, El, color=:blue, label="El", show = true)
end

function analyze(i)
    El0 = E0L[i]
    Elnew = UT.scale_for_inclusion_contact_point(El0, El) 
    plot_config!(El0, El, Elnew)
    println(El0∈El)

    Elnew = UT.scale_for_noninclusion_contact_point(El0, El) 
    plot_config!(El0, El, Elnew)
    println(UT.intersect(El0, El))
end

# We define some ellipsoids
c = [1.5; 1.5]
P = [4.0 0.5;       
     0.5 6.0]
El = UT.Ellipsoid(P, c)

P0 = [0.4 -0.1;
      -0.1 0.5]
vals = [4.1, 3.32, 2.8, 2.4]
E0L = [UT.Ellipsoid(P0, [c0x; c0x-0.2]) for c0x in vals]


# case 1 : non intersection
analyze(1)

# case 2 : non intersection
analyze(2)

# case 3 : intersection, non inclusion
analyze(3)

# case 4 : inclusion
analyze(4)




