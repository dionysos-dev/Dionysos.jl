module PWAsys

using Polyhedra
using FillArrays
using MathematicalSystems, HybridSystems
using StaticArrays
using LinearAlgebra
using Dionysos
using Dionysos.Control



function build_PWA_system(lib, dt, x_f, Usz)
    eye(n) = diagm(ones(n)) # I matrix
    # Define system
    N_region = 3
    n_sys = 2 
    n_u = 2; 
    #PWA partitions
    repX1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
    pX1 = polyhedron(repX1, lib);
    repX2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
    pX2 = polyhedron(repX2, lib);
    repX3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
    pX3 = polyhedron(repX3, lib);

    #control input bounded region
    if n_u>1 
        repU = intersect([HalfSpace(SVector{n_u}(-(1:n_u .==i)  ), Usz) ∩ HalfSpace(SVector{n_u}((1:n_u .==i)), Usz) for i in 1:n_u]...);
    else
        repU = HalfSpace(SVector{n_u}(-[1.0]), Usz) ∩ HalfSpace(SVector{n_u}([1.0]), Usz) 
    end
    pU = polyhedron(repU, lib);
    pX = [pX1 pX2 pX3];



    A = Vector{SMatrix{n_sys,n_sys,Float64}}(undef, N_region)
    B = Vector{SMatrix{(n_sys,n_u),Float64}}(undef, N_region)
    g = Vector{SVector{n_sys,Float64}}(undef, N_region)
    A[1] = SMatrix{2,2}(eye(n_sys)+[1 30 ;
        -10 1]*dt);
    A[2] = transpose(A[1]);
    A[3] = A[1];
    B = fill(SMatrix{2,2}(eye(n_u)*dt), N_region)
    g[1] = -SMatrix{2,1}([10; 10])*0.01;
    g[2] = g[1]*0;
    g[3] = -g[1];

    # automata for pwa switching between partitions (unecessary?)
    a = GraphAutomaton(N_region)
    add_transition!(a, 1, 2, 2); 
    add_transition!(a, 2, 1, 1);
    add_transition!(a, 1, 1, 1);
    add_transition!(a, 3, 2, 2);
    add_transition!(a, 2, 2, 2);
    add_transition!(a, 2, 3, 3);
    add_transition!(a, 3, 3, 3);

    # subsystems
    systems = [ConstrainedAffineControlDiscreteSystem(A[i], B[i], g[i], pX[i], pU) for i in 1:N_region];

    switching = AutonomousSwitching()
    switchings = fill(switching, 1)
    resetmaps = []
    
    system =  HybridSystem(a, systems, resetmaps, switchings)

    system.ext[:obstacles] = [Dionysos.Utils.HyperRectangle(SVector(0.0, -1.0), SVector(0.25, 1.5)), 
                              Dionysos.Utils.HyperRectangle(SVector(0.0, 1.25), SVector(1.0, 1.5))] 
    system.ext[:goal] = x_f; 
    return system # pwa system
end

""""
    build_problem(lib, dt=0.01, Usz=50, x_0 = [2.0,-2.0], x_f = [-2.0, 1.0], N = -1)

This function create the system with `PWAsys` and instantiates our OptimalControlProblem 
by defining the transition costs.
Notice that `state_cost` is defined to be zero for each mode/discrete state
of the system and the `transition_cost` is defined to be a quadratic function
of the state and the input.

Notice that we used `Fill` for all `N` time steps as we consider time-invariant costs.
"""
function build_problem(lib, dt=0.01, Usz=50, x_0 = [2.0,-2.0], x_f = [-2.0, 1.0], N = -1)
    system = build_PWA_system(lib, dt, x_f, Usz)
    n_sys = size(system.modes[1].A,1);
    n_u = size(system.modes[1].B,2);
    
    state_cost = ZeroFunction()
    transition_cost = Fill(QuadraticStateControlFunction(Matrix(I(n_sys)),Matrix(I(n_u)),zeros(n_sys,n_u),zeros(n_sys),zeros(n_u),0.0),nmodes(system))
    problem = OptimalControlProblem(
        system,
        Nothing, x_0,
        Fill(state_cost,nmodes(system)),
        Fill(transition_cost, ntransitions(system)),
        Nothing,
        N,
    )
    return problem
end

end