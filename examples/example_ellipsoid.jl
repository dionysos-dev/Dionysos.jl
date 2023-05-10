include("../src/Dionysos.jl")
using .Dionysos
UT = Dionysos.Utils
using Plots, Colors, LinearAlgebra

function check_inclusion(El0,El)
    println(El0∈El)

    Elnew = UT.scale_for_inclusion_contact_point(El0, El) 
    println(Elnew)

    p = plot(aspect_ratio=:equal)
    UT.plotE!(El0, color=:red, label="El0")
    UT.plotE!(Elnew, color=:green, label="Elnew")
    UT.plotE!(El, color=:blue, label="El")
    display(p)
end

function check_intersection(El0,El)
    println(UT.intersect(El0, El))

    Elnew = UT.scale_for_noninclusion_contact_point(El0, El) 
    println(Elnew)

    p = plot(aspect_ratio=:equal)
    UT.plotE!(El0, color=:red, label="El0")
    UT.plotE!(El, color=:blue, label="El")
    UT.plotE!(Elnew, color=:green, label="Elnew")
    display(p)
end

a = 2.5 #0.8 #1.2 #1.52 #2.5 
c0 = [1.6+a; 1.4+a]
P0 = [0.4 -0.1;
     -0.1 0.5]

c = [1.5; 1.5]
P = [4.0 0.5;       
     0.5 6.0]

El0 = UT.Ellipsoid(P0, c0)
El = UT.Ellipsoid(P, c)

# check_intersection(El0,El)
check_inclusion(El0,El)
