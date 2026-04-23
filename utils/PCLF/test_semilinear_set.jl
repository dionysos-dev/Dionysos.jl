using LazySets
using Dionysos
using Plots
const DI = Dionysos
const UT = DI.Utils

const Poly = LazySets.HPolytope

P1 = Poly([
    LazySets.HalfSpace([-1.144802508840479, -0.25374062921105384], 2.0867448918873683),
    LazySets.HalfSpace([0.5879672474030163, 0.4971867034225444], -0.7807957078098577),
    LazySets.HalfSpace([-0.624215261529317, -0.5910584452022201], 0.7806957078098578),
])

P2 = Poly([
    LazySets.HalfSpace([-0.491591697066812, -0.5068485705200367], 0.8984212500169478),
    LazySets.HalfSpace([0.5879672474030163, 0.4971867034225444], -0.7807957078098577),
    LazySets.HalfSpace([-0.5108516502549141, -0.43196828009757515], 0.6783964517549395),
])

println(LazySets.isempty(P1))
println(LazySets.isempty(P2))

#p = plot(; aspect_ratio = :equal);
#plot!(p,P1)
#plot!(p,P2)
#display(p)

R = Poly(vcat(LazySets.constraints_list(P1), LazySets.constraints_list(P2)))
println(LazySets.isempty(R))

R = UT.clean_poly(R)
