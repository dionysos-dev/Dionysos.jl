using Dionysos
using StaticArrays
using LinearAlgebra
using Plots

const DI = Dionysos
const UT = DI.Utils

function plot_config!(El0, El, Elnew)
    fig = plot(aspect_ratio=:equal)
    if Elnew∈El0
        plot!(fig, El0, color=:red, label="El0")
        plot!(fig, Elnew, color=:green, label="Elnew")
    else
        plot!(fig, Elnew, color=:green, label="Elnew")
        plot!(fig, El0, color=:red, label="El0")
    end
    plot!(fig, El, color=:blue, label="El")
    display(fig)
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

c = [1.5; 1.5]
P = [4.0 0.5;
     0.5 6.0]
El = UT.Ellipsoid(P, c)

P0 = [0.4 -0.1;
      -0.1 0.5]
vals = [4.1, 3.32, 2.8, 2.4]
E0L = [UT.Ellipsoid(P0, [c0x; c0x-0.2]) for c0x in vals]

analyze(1)

analyze(2)

analyze(3)

analyze(4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

