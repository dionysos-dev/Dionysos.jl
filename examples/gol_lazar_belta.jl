using FillArrays
using Polyhedra
using MathematicalSystems, HybridSystems
using SemialgebraicSets
using Dionysos

function gol_lazar_belta(lib)
    function rect(x_l, x_u)
        r = HalfSpace([-1], -x_l) ∩ HalfSpace([1], x_u)
        return polyhedron(r, lib)
    end
    function rect(x_l, x_u, y_l, y_u)
        r = HalfSpace([-1, 0], -x_l) ∩ HalfSpace([1, 0], x_u) ∩ HalfSpace([0, -1], -y_l) ∩ HalfSpace([0, 1], y_u)
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

    domains = [pC1, pC2, pC3, pC4, pC5, pC6, pC7, pC8, pT, pA, pB,
               pC1, pC2, pC3, pC4, pC5, pC6, pC7, pC8, pT]

    #Automaton:
    automaton = LightAutomaton(length(domains))

    #Resetmaps: (guards + reset)
    A = [1.0 1.0
         0.0 1.0]
    B = reshape([0.5, 1.0], 2, 1)

    reset_map = ConstrainedLinearControlMap(A, B, FullSpace(), pU)

    function back(from, to)
        # The `hrep` and `polyhedron` are workaround for an issue similar to
        # https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/216
        return polyhedron(hrep(project(polyhedron([A B] \ hrep(domains[to]) ∩ hrep(cU), lib), 1:2)) ∩ hrep(domains[from]), lib)
    end

    function is_enabled(from, to)
        all([0, 11]) do offset
            return !(from == offset + 3 && !((to - offset) in 1:3)) &&
            !(to == offset + 3 && !((from - offset) in 1:3))
        end
    end

    cU = polyhedron(HalfSpace([0, 0, -1], 2) ∩ HalfSpace([0, 0, 1], 2), lib)

    k = 0
    function maybe_add_transition(from, to)
        if is_enabled(from, to)
            if reachable_any_to_any(domains[from], domains[to], reset_map)
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
        Fill(reset_map, k),
        Fill(ControlledSwitching(), k)
    )

    system.ext[:q_A] = 10
    system.ext[:q_B] = 11
    system.ext[:q_T] = nmodes(system)
    system.ext[:obstacles] = [pO1, pO2]
    return system
end
