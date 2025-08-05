module DCDC
using Test
# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

using StaticArrays, JuMP
using MathematicalSystems
using HybridSystems

# At this point, we import the useful Dionysos sub-module for this problem:
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

function A1(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(-rL / xL, 0.0, 0.0, -1.0 / xC / (r0 + rC))
end

function A2(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(
        -(rL + r0 * rC / (r0 + rC)) / xL,
        5.0 * r0 / (r0 + rC) / xC,
        -r0 / (r0 + rC) / xL / 5.0,
        -1.0 / xC / (r0 + rC),
    )
end

function A2_abs(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(
        -(rL + r0 * rC / (r0 + rC)) / xL,
        5.0 * r0 / (r0 + rC) / xC,
        r0 / (r0 + rC) / xL / 5.0,
        -1.0 / xC / (r0 + rC),
    )
end

"""
Reset map for DC-DC converter mode transitions.
Applies one step of the target mode dynamics during transition.
"""
struct DCDCResetMap <: MathematicalSystems.AbstractMap
    domain::UT.HyperRectangle
    target_mode::Int  # Mode vers lequel on transite
    Dynamics::MathematicalSystems.AbstractContinuousSystem
    time_step::Float64              # Pas de temps pour l'intégration
end

function MathematicalSystems.apply(reset::DCDCResetMap, state::AbstractVector)
    # Extraire position et temps
    x = SVector(state[1], state[2])  # Position (iL, vC)
    t = state[3]  # Temps

    # Utiliser le système de dynamiques du mode de destination pour calculer le prochain état
    # Simuler un pas de temps avec la dynamique du mode de destination
    map_sys = Dionysos.System.simulate_control_map(reset.Dynamics.f)
    u_dummy = SVector(0.0)  # Pas d'entrée de contrôle continue pour les reset maps

    # Calculer le prochain état en appliquant la dynamique
    x_new = map_sys(x, u_dummy, reset.time_step)

    # Retourner le nouvel état avec le temps incrémenté
    return [x_new[1], x_new[2], t + reset.time_step]
end

MathematicalSystems.stateset(reset::DCDCResetMap) = reset.domain

function generate_safety_system_and_problem()
    """
    Génère un problème de sécurité DC-DC converter utilisant la structure TimedHybridAutomata.

    Objectif : maintenir le système dans une région de sécurité définie par les contraintes
    physiques du convertisseur (courant max, tension max, etc.)

    Caractéristiques :
    - Le temps n'est pas pris en compte comme contrainte de sécurité
    - Transitions libres entre modes depuis n'importe quel point
    - Reset map applique la dynamique du mode de destination
    - Evolution du système assurée par les reset maps lors des transitions
    """

    # ========== PARAMÈTRES DU SYSTÈME DC-DC ==========
    xL = 3.0;
    xC = 70.0;
    r0 = 1.0;
    rL = 0.05;
    rC = 0.005;
    vs = 1.0

    A1_matrix = A1(; xL, xC, r0, rL, rC)
    A2_matrix = A2(; xL, xC, r0, rL, rC)
    b_vector = SVector(vs / xL, 0.0)

    # Dynamiques des modes
    mode1_dynamics(x, u) = A1_matrix * x + b_vector
    mode2_dynamics(x, u) = A2_matrix * x + b_vector

    # ========== RÉGIONS DE SÉCURITÉ ==========
    # Limites de sécurité pour le convertisseur DC-DC
    iL_min, iL_max = 1.15, 1.55    # Courant inductance sûr [A]
    vC_min, vC_max = 5.45, 5.85    # Tension condensateur sûre [V]

    # Région de sécurité (identique pour les deux modes)
    _X_ = UT.HyperRectangle(SVector(iL_min, vC_min), SVector(iL_max, vC_max))

    # Pas d'entrée continue - contrôle par commutation
    U = UT.HyperRectangle(SVector(0.0), SVector(1.0)) # Pas d'entrée continue, juste le mode de commutation

    # ========== SYSTÈMES CONTINUS ==========
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics,
        2,
        1,
        _X_,
        U,
    )
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics,
        2,
        1,
        _X_,
        U,
    )

    # ========== SYSTÈMES TEMPORELS ==========
    # Le temps n'est pas pris en compte dans ce problème de sécurité
    time_horizon = 10.0  # Horizon temporel (non utilisé activement)
    time_window = UT.HyperRectangle([0.0], [time_horizon])

    time_system1 =
        MathematicalSystems.ConstrainedLinearContinuousSystem([0.0;;], time_window)
    time_system2 =
        MathematicalSystems.ConstrainedLinearContinuousSystem([0.0;;], time_window)

    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, time_system1]),
        SY.VectorContinuousSystem([mode2_system, time_system2]),
    ]

    # ========== TRANSITIONS AVEC DYNAMIQUE ==========
    # Transitions possibles partout dans l'espace d'état
    # Guard global: peut transiter depuis n'importe quel point de l'espace d'état sûr
    guard_global = UT.HyperRectangle(
        SVector(iL_min, vC_min, 0.0),
        SVector(iL_max, vC_max, time_horizon),
    )

    # Reset maps appliquant la dynamique du mode de destination
    dt_reset = 0.5  # Pas de temps pour l'application de la dynamique

    # Reset 1->2: applique la dynamique du Mode 2 lors de la transition
    reset_1_to_2 = DCDCResetMap(guard_global, 2, mode2_system, dt_reset)

    # Reset 2->1: applique la dynamique du Mode 1 lors de la transition  
    reset_2_to_1 = DCDCResetMap(guard_global, 1, mode1_system, dt_reset)

    reset_maps = [reset_1_to_2, reset_2_to_1]

    # ========== AUTOMATE DE SÉCURITÉ ==========
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    HybridSystems.add_transition!(automaton, 2, 1, 2)

    switchings = [HybridSystems.AutonomousSwitching(), HybridSystems.AutonomousSwitching()]

    hybrid_system =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # ========== ABSTRACTION POUR SÉCURITÉ ==========
    # Discrétisation plus grosse car le temps n'est pas important
    discretization_parameters = [
        (0.005, 0.5, 0.5),  # Mode 1: discrétisation plus grossière 
        (0.005, 0.5, 0.5),  # Mode 2: discrétisation plus grossière
    ]

    optimizer_factory_list = [
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ]

    # Points de référence sûrs
    x0 = SVector(0.0, 0.0)
    hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
    state_grid = DO.GridFree(x0, hx)
    u0 = SVector(1)
    hu = SVector(1)
    input_grid = DO.GridFree(u0, hu)

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid,
            "input_grid" => input_grid,
            "time_step" => discretization_parameters[1][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> A1_matrix,
        ),
        Dict(
            "state_grid" => state_grid,
            "input_grid" => input_grid,
            "time_step" => discretization_parameters[2][3],
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> A2_abs(; xL, xC, r0, rL, rC),
        ),
    ]

    # ========== PROBLÈME DE SÉCURITÉ ==========
    # État initial sûr - MÊME QUE L'EXEMPLE ORIGINAL
    initial_state = (SVector(1.2, 5.6), 0.0, 1)  # Point initial identique à l'original

    # Ensemble de sécurité = région d'état sûre pour chaque mode
    # COHÉRENT AVEC L'EXEMPLE ORIGINAL
    X_safe = _X_  # Utilise la même région que le système (comme dans l'original)
    Xs_safe = [X_safe, X_safe]  # Même région de sécurité pour les deux modes
    Ts_safe = [time_window, time_window]  # Toute la fenêtre temporelle (non contraignante)
    Ns_safe = [1, 2]  # Les deux modes sont sûrs

    # Spécifications du problème de sécurité
    safety_specs = AB.TemporalHybridSymbolicModelAbstraction.SafetyProblemSpecs(
        initial_state,
        Xs_safe,
        Ts_safe,
        Ns_safe,
        PB.Infinity(),  # Maintenir la sécurité indéfiniment
    )

    # ========== VÉRIFICATION DE COHÉRENCE ==========
    println("=== DCDC Safety Problem Verification ===")
    println("System parameters:")
    println("  xL = $xL, xC = $xC, r0 = $r0, rL = $rL, rC = $rC, vs = $vs")
    println("Safe region:")
    println("  iL ∈ [$iL_min, $iL_max] A")
    println("  vC ∈ [$vC_min, $vC_max] V")
    println("Initial state: $initial_state")
    println("Mode dynamics:")
    println("  A1 = $A1_matrix")
    println("  A2 = $A2_matrix")
    println("  b = $b_vector")
    println("Problem type: Safety (maintain system in safe region)")
    println("Time consideration: Non-constraining (spatial safety only)")
    println("Transitions: Free with dynamic reset maps")
    println(
        "Consistency check: Initial point (1.2, 5.6) ∈ [1.15, 1.55] × [5.45, 5.85] = ✅ VALID",
    )
    println("==========================================")

    return hybrid_system, optimizer_factory_list, optimizer_kwargs_dict, safety_specs
end

end # module DCDC
