module PWAsys

using StaticArrays, LinearAlgebra, Polyhedra
using MathematicalSystems, HybridSystems
using FillArrays, CDDLib

import Dionysos
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

function system(lib, dt, Usz, Wsz; simple = false)
    eye(n) = diagm(ones(n)) # I matrix
    # Define system
    N_region = 3
    n_sys = 2
    n_u = 2
    #PWA partitions
    repX1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
    pX1 = polyhedron(repX1, lib)
    repX2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
    pX2 = polyhedron(repX2, lib)
    repX3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
    pX3 = polyhedron(repX3, lib)

    #control input bounded region
    if n_u > 1
        repU = intersect(
            [
                HalfSpace(SVector{n_u}(-(1:n_u .== i)), Usz) ∩
                HalfSpace(SVector{n_u}((1:n_u .== i)), Usz) for i in 1:n_u
            ]...,
        )
    else
        repU = HalfSpace(SVector{n_u}(-[1.0]), Usz) ∩ HalfSpace(SVector{n_u}([1.0]), Usz)
    end
    pU = polyhedron(repU, lib)
    pX = [pX1 pX2 pX3]

    A = Vector{SMatrix{n_sys, n_sys, Float64}}(undef, N_region)
    B = Vector{SMatrix{(n_sys, n_u), Float64}}(undef, N_region)
    g = Vector{SVector{n_sys, Float64}}(undef, N_region)
    A[1] = SMatrix{2, 2}(eye(n_sys) + [
        1 30
        -10 1
    ] * dt)
    A[2] = transpose(A[1])
    A[3] = A[1]
    B = fill(SMatrix{2, 2}(eye(n_u) * dt), N_region)
    g[1] = -SMatrix{2, 1}([10; 10]) * 0.01
    g[2] = g[1] * 0
    g[3] = -g[1]

    # automata for pwa switching between partitions (unecessary?)
    a = GraphAutomaton(N_region)
    add_transition!(a, 1, 2, 2)
    add_transition!(a, 2, 1, 1)
    add_transition!(a, 1, 1, 1)
    add_transition!(a, 3, 2, 2)
    add_transition!(a, 2, 2, 2)
    add_transition!(a, 2, 3, 3)
    add_transition!(a, 3, 3, 3)

    # subsystems
    systems = [ConstrainedContinuousIdentitySystem(n_sys, pX[i]) for i in 1:N_region]

    switching = AutonomousSwitching()
    switchings = fill(switching, 1)
    resetmaps =
        [ConstrainedAffineControlMap(A[i], B[i], g[i], pX[i], pU) for i in 1:N_region]

    system = HybridSystem(a, systems, resetmaps, switchings)

    simple ? rectX = UT.HyperRectangle(SVector(-2.0, -1.5), SVector(-0.5, 1.3)) :
    rectX = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    simple ? obs = [] :
    obs = [
        UT.HyperRectangle(SVector(0.0, -1.0), SVector(0.25, 1.5)),
        UT.HyperRectangle(SVector(0.0, 1.25), SVector(1.0, 1.5)),
    ]
    Uaux = diagm(1:n_u)
    U = [(Uaux .== i) ./ Usz for i in 1:n_u] # matrices U_i
    W = Wsz * [
        -1 -1 1 1
        -1 1 -1 1
    ] * dt # polytope of disturbances 
    obs = UT.LazyUnionSetArray(obs)
    rectX = UT.LazySetMinus(rectX, obs)
    system.ext[:X] = rectX
    system.ext[:U] = U
    system.ext[:W] = W
    system.ext[:dt] = dt

    return system # pwa system
end

""""
    problem(lib, dt=0.01, Usz=50, x_0 = [2.0,-2.0], x_f = [-2.0, 1.0], N = -1)

This function create the system with `PWAsys` and instantiates our OptimalControlProblem 
by defining the transition costs.
Notice that `state_cost` is defined to be zero for each mode/discrete state
of the system and the `transition_cost` is defined to be a quadratic function
of the state and the input.

Notice that we used `Fill` for all `N` time steps as we consider time-invariant costs.

This problem was tackled in the paper [State-feedback Abstractions for Optimal Control of Piecewise-affine Systems](https://arxiv.org/abs/2204.00315).
"""
function problem(;
    lib = CDDLib.Library(),
    dt = 0.01,
    Usz = 50,
    Wsz = 5,
    x_0 = [1.8, -1.8],
    x_f = [-1.5, 1.0],
    N = Infinity(),
    simple = false,
)
    sys = system(lib, dt, Usz, Wsz; simple = simple)
    n_sys = size(sys.resetmaps[1].A, 1)
    n_u = size(sys.resetmaps[1].B, 2)
    if simple
        x_0 = [-1.0, -0.6]
        x_f = [-1.5, 1.0]
    end
    state_cost = UT.ZeroFunction()
    transition_cost = Fill(
        UT.QuadraticStateControlFunction(
            Matrix{Float64}(I(n_sys) * (dt^2)),
            Matrix{Float64}(I(n_u) * (dt^2)),
            zeros(n_sys, n_u),
            zeros(n_sys),
            zeros(n_u),
            0.0,
        ),
        nmodes(sys),
    )
    problem = OptimalControlProblem(
        sys,
        x_0,
        x_f,
        Fill(state_cost, nmodes(sys)),
        Fill(transition_cost, ntransitions(sys)),
        N,
    )
    return problem
end

end
