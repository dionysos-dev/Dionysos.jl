module GolLazarBelta

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem

using FillArrays
using Polyhedra
using MathematicalSystems, HybridSystems
using SemialgebraicSets

import CDDLib

# function to get vertices of polygon in the good cw order (for plotting)
function get_ordered_vertices(po)
    center = center_of_mass(po)
    p = [points(po)...]
    sort!(p; lt = (a, b) -> begin
        a_ang = atan((a - center)...) % 2pi
        b_ang = atan((b - center)...) % 2pi
        return a_ang < b_ang
    end)
    return p
end

function system(lib, T::Type)
    function rect(x_l, x_u)
        r = HalfSpace([-1], -T(x_l)) ∩ HalfSpace([1], T(x_u))
        return polyhedron(r, lib)
    end
    function rect(x_l, x_u, y_l, y_u)
        r =
            HalfSpace([-1, 0], -T(x_l)) ∩ HalfSpace([1, 0], T(x_u)) ∩
            HalfSpace([0, -1], -T(y_l)) ∩ HalfSpace([0, 1], T(y_u))
        return polyhedron(r, lib)
    end

    pX = rect(-10, 1.85, -10, 2)
    pU = rect(-2, 2)

    pT = rect(-0.5, 0.5, -0.5, 0.5)
    pO1 = rect(-10, -5, -10, -5)
    pO2 = rect(-5, 1.85, -4, -3)

    pA = rect(-6, -5, 1, 2)
    pB = rect(-5, -4, -3, -2)

    #Convex subsets:
    pC1 = rect(-6, -5, -5, 1)
    pC2 = rect(-10, -6, -5, 2)
    pC3 = rect(-5, 1.85, -10, -4)
    pC4 = rect(-4, -0.5, -3, 2)
    pC5 = rect(-5, -4, -2, 2)
    pC6 = rect(-0.5, 0.5, -3, -0.5)
    pC7 = rect(-0.5, 0.5, 0.5, 2)
    pC8 = rect(0.5, 1.85, -3, 2)

    domains = [
        pC1,
        pC2,
        pC3,
        pC4,
        pC5,
        pC6,
        pC7,
        pC8,
        pT,
        pA,
        pB,
        pC1,
        pC2,
        pC3,
        pC4,
        pC5,
        pC6,
        pC7,
        pC8,
        pT,
    ]

    #Automaton:
    automaton = GraphAutomaton(length(domains))

    #Resetmaps: (guards + reset)
    A = T[
        1 1
        0 1
    ]
    B = reshape(T[0.5, 1], 2, 1)

    #ResetMaps = Vector{ConstrainedLinearControlDiscreteSystem}(undef, length(domains)^2)
    cU = polyhedron(HalfSpace(T[0, 0, -1], 2) ∩ HalfSpace(T[0, 0, 1], 2), lib)

    function back(from, to)
        # The `hrep` and `polyhedron` are workaround for an issue similar to
        # https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/216
        return polyhedron(
            hrep(project(polyhedron([A B] \ hrep(domains[to]) ∩ hrep(cU), lib), 1:2)) ∩
            hrep(domains[from]),
            lib,
        )
    end

    function is_enabled(from, to)
        all([0, 11]) do offset
            return !(from == offset + 3 && !((to - offset) in 1:3)) &&
                   !(to == offset + 3 && !((from - offset) in 1:3))
        end
    end

    k = 0
    function maybe_add_transition(from, to)
        if is_enabled(from, to)
            guard = back(from, to)
            if !isempty(guard)
                k += 1
                add_transition!(automaton, from, to, k)
            end
        end
    end
    for from in 1:9, to in 1:11
        maybe_add_transition(from, to)
    end
    for from in 10:length(domains), to in 10:length(domains)
        maybe_add_transition(from, to)
    end

    system = HybridSystem(
        automaton,
        # Modes
        [ConstrainedContinuousIdentitySystem(2, i) for i in domains],
        # Reset maps
        Fill(ConstrainedLinearControlMap(A, B, FullSpace(), pU), k),
        Fill(ControlledSwitching(), k),
    )

    system.ext[:q_A] = 10
    system.ext[:q_B] = 11
    system.ext[:q_T] = nmodes(system)
    system.ext[:obstacles] = [pO1, pO2]
    return system
end

""""
    problem(lib, T::Type, q_0 = 3, x_0 = [1.0, -6.0], N = 11, zero_cost::Bool = true)

This function create the system with `GolLazarBelta.system`.

Then, we define initial conditions (continuous and discrete states) to this system
and set `N` as the search depth, i.e., the number of allowed time steps.

We instantiate our Optimal Control Problem by defining the state and transition costs.
Notice that `state_cost` is defined to be zero for each mode/discrete state
of the system and the `transition_cost` is defined to be `u_1^2` which is defined by the
quadratic form `u' * Q * u` with `Q = ones(1, 1)`.

Notice that we used `Fill` for all `N` time steps as we consider time-invariant costs.
"""
function problem(
    lib = CDDLib.Library(),
    T::Type = Float64;
    q_0 = 3,
    x_0 = [1.0, -6.0],
    N = 11,
    zero_cost::Bool = true,
)
    sys = system(lib, T)
    if zero_cost
        state_cost = Fill(UT.ZeroFunction(), nmodes(sys))
    else
        state_cost = [
            mode == sys.ext[:q_T] ? UT.ConstantFunction(zero(T)) :
            UT.ConstantFunction(one(T)) for mode in modes(sys)
        ]
    end
    transition_cost = UT.QuadraticControlFunction(ones(T, 1, 1))
    problem = PR.OptimalControlProblem(
        sys,
        (q_0, x_0),
        sys.ext[:q_T],
        Fill(state_cost, N),
        Fill(Fill(transition_cost, ntransitions(sys)), N),
        N,
    )
    return problem
end

end
