module TestNonLinearBemporadMorari

include("../GolLazarBelta/solvers.jl")
include("../../../problems/non_linear.jl")

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

using LinearAlgebra, Test
import CDDLib
import Polyhedra
using MathematicalSystems, HybridSystems, SemialgebraicSets
using FillArrays

function rect2d(lib, T, x_l, x_u, y_l, y_u)
    r =
        Polyhedra.HalfSpace([-1, 0], -T(x_l)) ∩ Polyhedra.HalfSpace([1, 0], T(x_u)) ∩
        Polyhedra.HalfSpace([0, -1], -T(y_l)) ∩ Polyhedra.HalfSpace([0, 1], T(y_u))
    return Polyhedra.polyhedron(r, lib)
end

"""
Build the MOI.ScalarNonlinearFunction expression for the non-linear dynamics:
    x1' = 1.1*x1 - 0.2*x2 - μ*x2^3 + Ts*u1
    x2' = 1.1*x2 + 0.2*x1 + μ*x1^3 + Ts*u2
"""
function nonlinear_dynamics(μ, Ts)
    snf(op, args...) = MOI.ScalarNonlinearFunction(op, Any[args...])
    return (x, u) -> [
        snf(
            :+,
            snf(:*, 1.1, x[1]),
            snf(:*, -0.2, x[2]),
            snf(:*, -μ, snf(:^, x[2], 3)),
            snf(:*, Ts, u[1]),
        ),
        snf(
            :+,
            snf(:*, 1.1, x[2]),
            snf(:*, 0.2, x[1]),
            snf(:*, μ, snf(:^, x[1], 3)),
            snf(:*, Ts, u[2]),
        ),
    ]
end

"""
Build a single-mode HybridSystem using the exact nonlinear dynamics
encoded as MOI.ScalarNonlinearFunction via NonlinearControlMap.
"""
function nonlinear_system(lib, T, μ, Ts, state_bounds, input_bounds)
    pX = rect2d(lib, T, state_bounds...)
    pU = rect2d(lib, T, input_bounds...)

    f = nonlinear_dynamics(μ, Ts)
    nlmap = ST.NonlinearControlMap(f, 2, 2)

    automaton = GraphAutomaton(1)
    add_transition!(automaton, 1, 1, 1)

    sys = HybridSystem(
        automaton,
        [ConstrainedContinuousIdentitySystem(2, pX)],
        Fill(nlmap, 1),
        Fill(ControlledSwitching(), 1),
    )
    return sys, pU
end

@testset "NonLinear BemporadMorari" begin
    T = Float64
    lib = CDDLib.Library()
    μ = T(0.00005)
    Ts = T(1.0)

    sys, pU =
        nonlinear_system(lib, T, μ, Ts, (-5.0, 5.0, -5.0, 5.0), (-10.0, 10.0, -10.0, 10.0))

    # Use Ipopt for NLP (nonlinear constraints)
    nlp_solver =
        optimizer_with_attributes(cont_solver.optimizer_constructor, MOI.Silent() => true)

    @testset "Depth $N" for (N, x_0) in [(1, T[-2.0, -2.0]), (3, T[-3.0, -3.0])]
        q_0 = 1
        q_T = 1

        state_cost = Fill(UT.ZeroFunction(), nmodes(sys))
        transition_cost = UT.QuadraticControlFunction(Matrix{T}(I, 2, 2))
        problem = PR.OptimalControlProblem(
            sys,
            (q_0, x_0),
            q_T,
            Fill(state_cost, N),
            Fill(Fill(transition_cost, ntransitions(sys)), N),
            N,
        )

        algo = optimizer_with_attributes(
            OP.BemporadMorari.Optimizer{T},
            "continuous_solver" => nlp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0,
        )

        optimizer = MOI.instantiate(algo)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
        MOI.optimize!(optimizer)

        @test MOI.get(optimizer, MOI.TerminationStatus()) in
              [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]

        xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute())
        @test length(xu.x) == N
        @test length(xu.u) == N

        # Verify dynamics: x_next ≈ f(x_prev, u)
        for t in 1:N
            x_prev = t == 1 ? x_0 : xu.x[t - 1]
            x_next = xu.x[t]
            u_t = xu.u[t]
            f_val = [
                1.1 * x_prev[1] - 0.2 * x_prev[2] - μ * x_prev[2]^3 + Ts * u_t[1],
                1.1 * x_prev[2] + 0.2 * x_prev[1] + μ * x_prev[1]^3 + Ts * u_t[2],
            ]
            @test x_next ≈ f_val atol = 1e-4
        end

        obj = MOI.get(optimizer, MOI.ObjectiveValue())
        @test isfinite(obj)
        @test obj >= -1e-6
    end
end

end # module
