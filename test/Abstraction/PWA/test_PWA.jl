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


# @testset "PWA_LMIs" begin
#     eye(n) = diagm(ones(n))
#     n = 2
#     m = 1
#     A = [0.7    0.6
#         0.5    0.4];

#     B = hcat([.2; .5]);

#     H = 0.25*eye(n);

#     g = [0.2; 0.3];

#     c = [1;1];

#     cp = [3;5];


#     P = eye(n);
#     Pp = P;


#     W = H*[-1 -1  1 1;
#          -1  1 -1 1];

#     L = [eye(n+m) zeros(n+m,1)];
#     U = [1.0/(10)];

#     repX = HalfSpace(SVector{2}([-1.0, 0.0]), 6) ∩  HalfSpace(SVector{2}([1.0, 0.0]), 6) ∩ #
#            HalfSpace(SVector{2}([0.0, -1.0]), 6) ∩  HalfSpace(SVector{2}([0.0, 1.0]), 6);
            
#     pX = polyhedron(repX, lib);
#     repU = HalfSpace(SVector{1}([-1.0]), inv(U[1])) ∩ HalfSpace(SVector{1}([1.0]), inv(U[1]));
#     pU = polyhedron(repU, lib);

#     subsys = ConstrainedAffineControlDiscreteSystem(A, B, g, pX, pU)
    
#     res, cost, kappa = AB._has_transition(subsys,P,c,Pp,cp,W,L,U)
#     #println(res, cost, kappa)


#     @test res == true
#     @test cost ≈ 75.50 atol=1e-2


#     res, cost, kappa = AB._has_transition(subsys,P,c,Pp,cp*2,W,L,U)

#     @test res == false
    
# end




@testset "PWA_Trans" begin
    eye(n) = diagm(ones(n))
    # Define system
    N_region = 3
    n_sys = 2 
    n_u = 2; 
    Uaux = diagm(1:n_u)
    Usz = 1
    U = [(Uaux.==i)./Usz for i in 1:n_u];

    repX1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
    pX1 = polyhedron(repX1, lib);
    repX2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
    pX2 = polyhedron(repX2, lib);
    repX3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
    pX3 = polyhedron(repX3, lib);
    #repU = HalfSpace(SVector{1}([-1.0]), 1/U[1]) ∩ HalfSpace(SVector{1}([1.0]), 1/U[1]);
    repU = intersect([HalfSpace(SVector{n_u}(-(1:n_u .==i)  ), Usz) ∩ HalfSpace(SVector{n_u}((1:n_u .==i)), Usz) for i in 1:n_u]...);

    pU = polyhedron(repU, lib);
    pX = [pX1 pX2 pX3];



    W = [-1 -1  1 1;
         -1  1 -1 1]*0.001;

    PWAdomains = [pX1, pX2, pX3];

    A = Vector{SMatrix{n_sys,n_sys,Float64}}(undef, N_region)
    B = Vector{SMatrix{(n_sys,1),Float64}}(undef, N_region)
    H = Vector{SMatrix{(n_sys,2),Float64}}(undef, N_region)
    g = Vector{SVector{n_sys,Float64}}(undef, N_region)
    A[1] = SMatrix{2,2}([0.2 0.6 ;
        -0.6 0.2]*0.1);
    A[2] = -A[1];
    A[3] = 0.5*A[1];
#    B = fill(SMatrix{2,1}([0.3; 0.4]*10), N_region)
    B = fill(SMatrix{2,2}(eye(n_u)), N_region)
    H = fill(SMatrix{2,2}([0.1 0.0; 0.0 0.1]), N_region)
    g[1] = SMatrix{2,1}([0.0; 1.0]);
    g[2] = g[1]*0;
    g[3] = -g[1];

    a = LightAutomaton(N_region)
    add_transition!(a, 1, 2, 2); 
    add_transition!(a, 2, 1, 1);
    add_transition!(a, 1, 1, 1);
    add_transition!(a, 3, 2, 2);
    add_transition!(a, 2, 2, 2);
    add_transition!(a, 2, 3, 3);
    add_transition!(a, 3, 3, 3);

    systems = [ConstrainedAffineControlDiscreteSystem(A[i], B[i], g[i], pX[i], pU) for i in 1:N_region];



    switching = AutonomousSwitching()
    switchings = fill(switching, 1)


    resetmaps = []


    system = HybridSystem(a, systems, resetmaps, switchings)
    # Create Abstraction



    rectX = AB.HyperRectangle(SVector(-2.0, -2), SVector(2.0, 2));
    rectU = AB.HyperRectangle(SVector(-2.0, -2, -2), SVector(2.0, 2, 2));

    X_origin = SVector(0.0, 0.0);
    X_step = SVector(1.0/3, 1.0/3);
    Xgrid = AB.GridRectangular(X_origin, X_step, rectX);

    domainX = AB.DomainList(Xgrid);
    AB.add_set!(domainX, rectX, AB.OUTER)


    U_origin = SVector(0.0, 0.0, 0.0);
    U_step = SVector(1.0/5, 1.0/5, 10/5);
    Ugrid = AB.GridFree(U_origin, U_step);

    domainU = AB.DomainList(Ugrid);
    AB.add_set!(domainU, rectU, AB.INNER)


    symmodel = AB.NewSymbolicModelListList(domainX, domainU);



    L = [eye(n_sys+n_u) zeros(n_sys+n_u,1)];

    transitionCost = Dict()
    transitionKappa = Dict()
    AB.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionKappa, system, W, L, U)

    
#    @test true == false
    
end

sleep(0.1) # used for good printing
println("End test")

end # module TestMain