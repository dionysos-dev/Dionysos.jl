#include("../../../src/Abstraction/PWA/ellipsoidal_transitions.jl")
include("../../../src/Abstraction/abstraction.jl")

module TestMain

using Test
using LinearAlgebra
using Polyhedra
using StaticArrays
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets

lib = CDDLib.Library()
using Main.Abstraction
AB = Main.Abstraction


sleep(0.1) # used for good printing
println("Started test")


@testset "PWA_LMIs" begin
    eye(n) = diagm(ones(n))
    n = 2
    m = 1
    A = [0.7    0.6
        0.5    0.4];

    B = [.2; .5];

    H = 0.25*eye(n);

    g = [0.2; 0.3];

    c = [1;1];

    cp = [3;5];


    P = eye(n);
    Pp = P;


    W = [-1 -1  1 1;
         -1  1 -1 1];

    L = [eye(n+m) zeros(n+m,1)];
    U = [1.0/(10)];

    repX = HalfSpace(SVector{2}([-1.0, 0.0]), 6) ∩  HalfSpace(SVector{2}([1.0, 0.0]), 6) ∩ #
           HalfSpace(SVector{2}([0.0, -1.0]), 6) ∩  HalfSpace(SVector{2}([0.0, 1.0]), 6);
            
    pX = polyhedron(repX, lib);
    repU = HalfSpace(SVector{3}([-1.0,0,0]), 10) ∩ HalfSpace(SVector{3}([1.0,0,0]), 10);
    pU = polyhedron(repU, lib);

    subsys = ConstrainedAffineControlDiscreteSystem(A, [B H], g, pX, pU)
    
    res, cost, kappa = AB._has_transition(subsys,P,c,Pp,cp,W,L,U)
    #println(res, cost, kappa)


    @test res == true
    @test cost ≈ 75.50 atol=1e-2


    res, cost, kappa = AB._has_transition(subsys,P,c,Pp,cp*2,W,L,U)

    @test res == false
    
end

sleep(0.1) # used for good printing
println("End test")

end # module TestMain