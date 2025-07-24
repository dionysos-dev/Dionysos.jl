module TestMain

using Test
using StaticArrays, Plots, MathematicalSystems, HybridSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
  
sleep(0.1) 
println("Started test")

@testset "TimedHybridAutomata - simple, Time taken into account" begin
    #
    # Test 1 : Système hybride à 2 modes (état 1D + temps 1D)
    # - Mode 1 : dx/dt = 0.5*x + u, x ∈ [-1,1], u ∈ [-0.5,0.5], temps ∈ [0,2]
    # - Mode 2 : dx/dt = 0.8*x + u, x ∈ [-2,2], u ∈ [-0.3,0.3], temps ∈ [3,5]
    # - 1 transition (mode 1 → mode 2) avec une reset map qui force x = 0.0, t = 3.0
    # - Reset map : FixedPointResetMap (ignore l'état d'entrée, renvoie toujours la cible)  
    #
    # 

    # === MODE 1 ===
    # Système Black Box pour le mode 1 : dx/dt = 0.5*x + u
    function mode1_dynamics(x, u)
        return [0.5 * x[1] + u[1]]
    end
    # Domaines pour le mode 1
    X1 = UT.HyperRectangle([-1.0], [1.0])   # État x ∈ [-1,1]
    U1 = UT.HyperRectangle([-0.5], [0.5])   # Entrée u ∈ [-0.5,0.5]
                
    # Système Black Box du mode 1
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics, 1, 1, X1, U1
    )

    # === MODE 2 ===
    # Système Black Box pour le mode 2 : dx/dt = 0.8*x + u
    function mode2_dynamics(x, u)
        return [0.8 * x[1] + u[1]]
    end

    # Domaines pour le mode 2
    X2 = UT.HyperRectangle([-2.0], [2.0])   # État x ∈ [-2,2]
    U2 = UT.HyperRectangle([-0.3], [0.3])   # Entrée u ∈ [-0.3,0.3]

    # Système Black Box du mode 2
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics, 1, 1, X2, U2
    )
    # === DÉFINITION TEMPS DES MODES ===
    # Création du modèle temporel pour le mode 1
    temps_mode1 = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([0.0], [2.0]))

    # Création du modèle temporel pour le mode 2
    temps_mode2 = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([3.0], [5.0]))

    # === RESET MAP ===
    # Reset map simple : x2 = 0.0 quand on passe du mode 1 au mode 2
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end

    function MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector)
        # On ignore l'état d'entrée, on renvoie toujours la cible
        return reset.target
    end

    function MathematicalSystems.stateset(reset::FixedPointResetMap)
        # Le domaine de la reset map (ici, l'hyperrectangle du mode 1)
        return reset.domain
    end

    guard_1 = UT.HyperRectangle([0.0,1.0], [1.0 ,2.0])  # Condition de garde pour le reset
    reset_map = FixedPointResetMap(guard_1, [0.0, 3.0])

    # === AUTOMATE ===
    # Création de l'automate avec 2 états
    automaton = HybridSystems.GraphAutomaton(2)

    # Transition du mode 1 vers le mode 2
    HybridSystems.add_transition!(automaton, 1, 2, 1)

    # === CONSTRUCTION DU SYSTÈME HYBRIDE ===
    # Systèmes par mode
    modes_systems = [SY.VectorContinuousSystem([mode1_system, temps_mode1]), SY.VectorContinuousSystem([mode2_system, temps_mode2])]

    # Reset maps (un pour chaque transition)
    reset_maps = [reset_map]

    # Mécanismes de switching (autonomous pour simplicité)
    switchings = [HybridSystems.AutonomousSwitching()]

    # Construction du système hybride
    hs = HybridSystems.HybridSystem(
        automaton,
        modes_systems,
        reset_maps,
        switchings
    )

    # Bornes de croissance (matrices 1x1 ici, mais peut être plus grand si multidim)
    growth_bounds = SVector(SMatrix{1,1}(0.5), SMatrix{1,1}(0.8))

    # Paramètres de discrétisation pour chaque mode
    param_discretisation = [
        (0.5, 0.5, 0.5),  # (dx, du, dt) pour le mode 1
        (0.5, 0.3, 0.5)   # (dx, du, dt) pour le mode 2
    ]
    # Construction du modèle symbolique hybride temporel (prise en compte du temps)
    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,  # Matrices de croissance pour chaque mode
        param_discretisation  # Paramètres de discrétisation (état, entrée, temps) pour chaque mode
    )

    # ===================== TESTS DE ROBUSTESSE =====================
    
    # GlobalInputMap cohérent
    gim = hybrid_model.global_input_map
    @test gim.switching_inputs == 1
    @test length(gim.continuous_to_global) == 6   
    @test length(gim.switching_to_global) == 1
    # Fonctions d'accès GlobalInputMap
    @test SY.TimedHybridAutomata.get_global_input_id(gim, 1, 1) == 1
    @test SY.TimedHybridAutomata.get_global_input_id(gim, 2, 1) == 4
    @test SY.TimedHybridAutomata.get_switching_global_id(gim, 1) == 7
    @test SY.TimedHybridAutomata.is_continuous_input(gim, 1) == true
    @test SY.TimedHybridAutomata.is_switching_input(gim, gim.switching_range.start) == true
    # get_local_input_info
    typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, 1)
    @test typ == :continuous
    typ2, info2 = SY.TimedHybridAutomata.get_local_input_info(gim, gim.switching_range.start)
    @test typ2 == :switching

    # Modèle symbolique dynamique : nombre d'états > 0
    for sm in hybrid_model.symmodels
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(sm) == 3
    end
    # Modèle symbolique de temps : pas de temps correct
    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end


    autom = hybrid_model.autom
    #  Les états de l'automate sont bien dans int2aug_state
    for i in 1:length(hybrid_model.int2aug_state)
        @test haskey(hybrid_model.aug_state2int, hybrid_model.int2aug_state[i])
    end


    sm = hybrid_model.symmodels[1]
    tm = hybrid_model.time_symbolic_models[1]
    # find_symbolic_state : état existant
    x0 = [0.0]
    idx = SY.TimedHybridAutomata.find_symbolic_state(sm, x0)
    @test idx > 0
    # find_time_index : temps existant
    t0 = tm.tsteps[1]
    tidx = SY.TimedHybridAutomata.find_time_index(tm, t0)
    @test tidx == 1
    # extract_spatial_part et extract_temporal_part
    guard = guard_1
    sp = SY.TimedHybridAutomata.extract_spatial_part(guard)
    tp = SY.TimedHybridAutomata.extract_temporal_part(guard)
    @test sp.lb == guard.lb[1:end-1]
    @test tp == [guard.lb[end], guard.ub[end]]
    # get_time_indices_from_interval
    indices = SY.TimedHybridAutomata.get_time_indices_from_interval(tm, tp)
    @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    # Cas d'erreur : find_symbolic_state sur état hors domaine
    idx_bad = SY.TimedHybridAutomata.find_symbolic_state(sm, [100.0])
    @test idx_bad == 0
    
end

@testset "TimedHybridAutomata - simple, Time not taken into account" begin
    #
    # Test 2 : Système hybride à 2 modes (état 1D + temps 1D, temps figé)
    # - Mode 1 : dx/dt = 0.5*x + u, x ∈ [-1,1], u ∈ [-0.5,0.5], temps ∈ [0,2]
    # - Mode 2 : dx/dt = 0.8*x + u, x ∈ [-2,2], u ∈ [-0.3,0.3], temps ∈ [3,5]
    # - 1 transition (mode 1 → mode 2) avec une reset map qui force x = 0.0, t = 0.0
    # - Reset map : FixedPointResetMap (ignore l'état d'entrée, renvoie toujours la cible)
    #
    # === MODE 1 ===
    # Définition du système du mode 1 : dynamique simple dx/dt = 0.5*x + u
    function mode1_dynamics(x, u)
        return [0.5 * x[1] + u[1]]
    end
    # Domaine d'état x ∈ [-1,1] et domaine d'entrée u ∈ [-0.5,0.5] pour le mode 1
    X1 = UT.HyperRectangle([-1.0], [1.0])
    U1 = UT.HyperRectangle([-0.5], [0.5])
    # Système Black Box du mode 1
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics, 1, 1, X1, U1
    )

    # === MODE 2 ===
    # Définition du système du mode 2 : dynamique dx/dt = 0.8*x + u
    function mode2_dynamics(x, u)
        return [0.8 * x[1] + u[1]]
    end
    # Domaine d'état x ∈ [-2,2] et domaine d'entrée u ∈ [-0.3,0.3] pour le mode 2
    X2 = UT.HyperRectangle([-2.0], [2.0])
    U2 = UT.HyperRectangle([-0.3], [0.3])
    # Système Black Box du mode 2
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics, 1, 1, X2, U2
    )

    # === DÉFINITION DES SYSTÈMES DE TEMPS ===
    # Modèle temporel pour chaque mode (ici, dynamique figée à 0)
    temps_mode1 = MathematicalSystems.ConstrainedLinearContinuousSystem([0.0;;], UT.HyperRectangle([0.0], [2.0]))
    temps_mode2 = MathematicalSystems.ConstrainedLinearContinuousSystem([0.0;;], UT.HyperRectangle([3.0], [5.0]))

    # === RESET MAP ===
    # Définition d'une reset map qui force l'état à [0.0, 0.0] lors du passage du mode 1 au mode 2
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end

    function MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector)
        # Ignore l'état d'entrée, renvoie toujours la cible
        return reset.target
    end

    function MathematicalSystems.stateset(reset::FixedPointResetMap)
        # Domaine de validité de la reset map
        return reset.domain
    end

    guard_1 = UT.HyperRectangle([0.0,-1.0], [1.0 ,1.0])  # Garde : x ∈ [0,1], t ∈ [-1,1] par défaut 
    reset_map = FixedPointResetMap(guard_1, [0.0, 0.0])

    # === AUTOMATE ===
    # Création d'un automate à 2 états (modes)
    automaton = HybridSystems.GraphAutomaton(2)
    # Ajout d'une transition du mode 1 vers le mode 2
    HybridSystems.add_transition!(automaton, 1, 2, 1)

    # === CONSTRUCTION DU SYSTÈME HYBRIDE ===
    # Construction des systèmes par mode (état + temps)
    modes_systems = [SY.VectorContinuousSystem([mode1_system, temps_mode1]), SY.VectorContinuousSystem([mode2_system, temps_mode2])]
    # Liste des reset maps (une par transition)
    reset_maps = [reset_map]
    # Mécanismes de switching (ici, switching autonome)
    switchings = [HybridSystems.AutonomousSwitching()]

    # Construction du système hybride (2 modes, 1 transition, pas de temps pris en compte)
    hs = HybridSystems.HybridSystem(
        automaton,
        modes_systems,
        reset_maps,
        switchings
    )

    # Bornes de croissance (matrices 1x1 pour chaque mode)
    growth_bounds = SVector(SMatrix{1,1}(0.5), SMatrix{1,1}(0.8))

    # Paramètres de discrétisation (état, entrée, temps) pour chaque mode
    param_discretisation = [
        (0.5, 0.5, 0.5),  # mode 1
        (0.5, 0.3, 0.5)   # mode 2
    ]

    # Construction du modèle symbolique hybride temporel (sans prise en compte du temps)
    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretisation
    )

    # ===================== TESTS DE ROBUSTESSE =====================
    # Vérifie que le modèle hybride a été construit correctement avec le temps pas pris en compte

    # Vérifie que tous les états augmentés ont time_index == 1
    for aug_state in hybrid_model.int2aug_state
        # aug_state = (state_symbol, time_symbol, mode_id)
        @test aug_state[2] == 1
    end

    # Vérifie que chaque modèle symbolique de temps n'a qu'un seul pas de temps
    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) == 1
        @test tm.is_active == false
    end

    
    # Vérifie que find_time_index retourne toujours 1 pour tout temps dans le domaine
    for tm in hybrid_model.time_symbolic_models
        tmin, tmax = tm.tsteps[1], tm.tsteps[1]
        for tval in [tmin, tmax, (tmin + tmax)/2]
            @test SY.TimedHybridAutomata.find_time_index(tm, tval) == 1
        end
    end

    # Vérifie que get_time_indices_from_interval retourne [1] pour tout intervalle valide
    for tm in hybrid_model.time_symbolic_models
        tmin, tmax = tm.tsteps[1], tm.tsteps[1]
        @test SY.TimedHybridAutomata.get_time_indices_from_interval(tm, [tmin, tmax]) == [1]
    end
end

@testset "TimedHybridAutomata - 3 modes, transitions complètes, reset maps variées" begin
    #
    # Test 3 : Système hybride à 3 modes (état 2D + temps 1D)
    # - Mode A : dx/dt = A_A*x + B_A*u, x ∈ [-2,2]^2, u ∈ [-1,1]^2, temps ∈ [0,3]
    # - Mode B : dx/dt = A_B*x + B_B*u, x ∈ [-2,2]^2, u ∈ [-2,1]^2, temps ∈ [1.5,2.5]
    # - Mode C : dx/dt = A_C*x + B_C*u, x ∈ [-2,2]^2, u ∈ [-1,1]^2, temps ∈ [2,3]
    # - 6 transitions (toutes directions) avec reset maps variées :
    #   * offset (translation), perm (permutation des composantes), tmin (borne inf. sur t)
    #   * chaque reset map borne x_new par la borne sup. du mode cible
    # - Reset map : DenseResetMap (permutation, offset, tmin, x_max)
    #
    # === Dynamiques linéaires 2D (matrices) ===
    A_A = [0.5 0.0; 0.0 0.6]
    B_A = [1.0 0.0; 0.0 1.0]
    A_B = [0.4 0.1; 0.1 0.7]
    B_B = [1.0 0.0; 0.0 1.0]
    A_C = [0.8 0.1; 0.1 0.9]
    B_C = [1.0 0.0; 0.2 1.0]

    # === Domaines des états et entrées ===
    X_A = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_A = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])
    X_B = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_B = UT.HyperRectangle([-2.0, -2.0], [1.0, 1.0])
    X_C = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_C = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])

    # === Systèmes continus par mode (black box linéaire) ===
    function dynA(x, u)
        return A_A * x + B_A * u
    end
    function dynB(x, u)
        return A_B * x + B_B * u
    end
    function dynC(x, u)
        return A_C * x + B_C * u
    end
    modeA = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynA, 2, 2, X_A, U_A)
    modeB = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynB, 2, 2, X_B, U_B)
    modeC = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynC, 2, 2, X_C, U_C)

    # === Systèmes temporels (1D, même pour tous) ===
    temps_A = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([0.0], [3.0]))
    temps_B = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([1.5], [2.5]))
    temps_C = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([2.0], [3.0]))

    # === VectorContinuousSystem par mode ===
    vcA = SY.VectorContinuousSystem([modeA, temps_A])
    vcB = SY.VectorContinuousSystem([modeB, temps_B])
    vcC = SY.VectorContinuousSystem([modeC, temps_C])

    # === Automate à 3 états ===
    automaton = HybridSystems.GraphAutomaton(3)
    # 6 transitions (toutes directions)
    HybridSystems.add_transition!(automaton, 1, 2, 1) # A→B
    HybridSystems.add_transition!(automaton, 2, 1, 2) # B→A
    HybridSystems.add_transition!(automaton, 1, 3, 3) # A→C
    HybridSystems.add_transition!(automaton, 3, 1, 4) # C→A
    HybridSystems.add_transition!(automaton, 2, 3, 5) # B→C
    HybridSystems.add_transition!(automaton, 3, 2, 6) # C→B


    # === Reset maps variées ===
    # Une reset map définit comment l'état change lors d'une transition entre deux modes.
    # Ici, DenseResetMap permet de modéliser des resets complexes :
    #  - domain : domaine de validité de la reset (garde associée à la transition)
    #  - offset : translation appliquée à la partie état (x)
    #  - perm   : permutation des composantes de l'état (ex : [2,1] échange x₁ et x₂)
    #  - tmin   : valeur minimale imposée à la composante temps après reset
    struct DenseResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle  # Domaine sur lequel la reset est valide (garde)
        offset::Vector{Float64}    # Décalage appliqué à l'état (x)
        perm::Vector{Int}          # Permutation des indices de l'état (x)
        tmin::Float64              # Temps minimal après reset
        x_max::Vector{Float64}     # Valeur max autorisée pour chaque composante de x après reset (borne supérieure du mode cible)
    end

    # Fonction d'application de la reset map :
    # - Prend l'état courant (x₁, x₂, t)
    # - Applique la permutation perm à x
    # - Ajoute l'offset à x permuté
    # - Met à jour t en imposant t >= tmin
    # - Retourne le nouvel état [x_new; t_new]
    function MathematicalSystems.apply(reset::DenseResetMap, state::AbstractVector)
        x = state[1:2]  # Partie état (supposée dimension 2)
        t = max(state[3], reset.tmin)  # Temps borné par tmin
        x_perm = x[reset.perm]         # Permutation des composantes
        # Application de l'offset puis bornage par x_max (composante par composante)
        x_new = min.(x_perm .+ reset.offset, reset.x_max)
        return vcat(x_new, t)          # Nouvel état concaténé
    end

    # Domaine de validité de la reset map (utilisé pour la vérification des transitions)
    function MathematicalSystems.stateset(reset::DenseResetMap)
        return reset.domain
    end


    # Guards (domaines de déclenchement des transitions) :
    # Chaque garde est un hyperrectangle sur (x₁, x₂, t) qui définit quand la transition est possible.
    # Les domaines sont cohérents avec les domaines d'états et temps définis plus haut pour chaque mode.
    guards = [
        UT.HyperRectangle([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),   # A→B : x₁∈[1,2], x₂∈[0,1], t∈[0,1]
        UT.HyperRectangle([0.0, 1.0, 1.5], [1.0, 2.0, 2.5]),   # B→A : x₁∈[0,1], x₂∈[1,2], t∈[1.5,2.5]
        UT.HyperRectangle([1.0, -1.0, 0.0], [2.0, 0.0, 1.0]),  # A→C : x₁∈[1,2], x₂∈[-1,0], t∈[0,1]
        UT.HyperRectangle([-1.0, 1.0, 2.0], [0.0, 2.0, 3.0]),  # C→A : x₁∈[-1,0], x₂∈[1,2], t∈[2,3]
        UT.HyperRectangle([1.0, 1.0, 1.5], [2.0, 2.0, 2.5]),   # B→C : x₁∈[1,2], x₂∈[1,2], t∈[1.5,2.5]
        UT.HyperRectangle([-1.0, -1.0, 2.0], [0.0, 0.0, 3.0])  # C→B : x₁∈[-1,0], x₂∈[-1,0], t∈[2,3]
    ]

    # Offsets et permutations pour chaque reset map :
    # - offsets : translation appliquée à (x₁, x₂) lors du reset
    # - perms   : permutation des composantes (ex : [2,1] échange x₁ et x₂)
    # - tmins   : valeur minimale imposée à t après reset
    #
    # ⚠️ Pour garantir la cohérence, il faut que pour tout x dans la garde,
    #    x[perm] + offset ∈ domaine d'état du mode cible (ex : X_B pour A→B).
    #    Sinon, le reset peut produire un état hors du domaine autorisé.
    offsets = [
        [0.1, 0.0],   # A→B : translation
        [0.0, -0.1],  # B→A : translation négative
        [0.0, 0.0],   # A→C : identité
        [0.0, 0.2],  # C→A : translation mixte
        [0.0, 0.0],   # B→C : identité
        [0.3, 0.0]    # C→B : translation positive
    ]
    perms = [
        [1,2],  # A→B : identité
        [2,1],  # B→A : permutation
        [1,2],  # A→C : identité
        [2,1],  # C→A : permutation
        [1,2],  # B→C : identité
        [2,1]   # C→B : permutation
    ]
    tmins = [1.5, 0.0, 2.0, 0.0, 2.0, 1.5]

    # Chaque reset map est associée à une transition et définit comment l'état change lors du saut.
    # Par exemple, la transition B→A applique une permutation et une translation négative, et impose t >= 0.
    #
    # Pour garantir la sûreté, il est conseillé de vérifier que l'offset ne fait pas sortir l'état du domaine cible.
    # Voir la fonction utilitaire ci-dessous.

    # Fonction utilitaire pour extraire la borne supérieure (ub) d'un HyperRectangle

    x_maxs = [X_B.ub, X_A.ub, X_C.ub, X_A.ub, X_C.ub, X_B.ub]  # transitions : 1→2, 2→1, 1→3, 3→1, 2→3, 3→2
    reset_maps = [
        DenseResetMap(guards[1], offsets[1], perms[1], tmins[1], x_maxs[1]),
        DenseResetMap(guards[2], offsets[2], perms[2], tmins[2], x_maxs[2]),
        DenseResetMap(guards[3], offsets[3], perms[3], tmins[3], x_maxs[3]),
        DenseResetMap(guards[4], offsets[4], perms[4], tmins[4], x_maxs[4]),
        DenseResetMap(guards[5], offsets[5], perms[5], tmins[5], x_maxs[5]),
        DenseResetMap(guards[6], offsets[6], perms[6], tmins[6], x_maxs[6])
    ]
    switchings = [HybridSystems.AutonomousSwitching() for _ in 1:6]

    modes_systems = [vcA, vcB, vcC]

    # === Construction du système hybride ===
    hs = HybridSystems.HybridSystem(
        automaton,
        modes_systems,
        reset_maps,
        switchings
    )

    # === Paramètres pour abstraction symbolique ===
    growth_bounds = SVector(
        SMatrix{2,2}(abs.(A_A)),
        SMatrix{2,2}(abs.(A_B)),
        SMatrix{2,2}(abs.(A_C))
    )
    param_discretisation = [
        (0.5, 0.5, 0.5),  # mode A
        (0.5, 0.5, 0.5),  # mode B
        (0.5, 0.5, 0.5)   # mode C
    ]

    # === Construction du modèle symbolique hybride temporel ===
    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretisation
    )
    # ===================== TESTS DE ROBUSTESSE =====================
    # Vérification de la cohérence des entrées globales
    # ===================== TESTS DE ROBUSTESSE =====================

    # 1. Vérification de la structure du modèle
    @test length(hybrid_model.symmodels) == 3
    @test length(hybrid_model.time_symbolic_models) == 3
    @test length(hybrid_model.int2aug_state) == length(hybrid_model.aug_state2int)

    # 2. Vérification du GlobalInputMap
    gim = hybrid_model.global_input_map
    @test gim.switching_inputs == 6
    @test length(gim.continuous_to_global) == 25 + 25 + 7*7  
    @test length(gim.switching_to_global) == 6
    # Fonctions d'accès GlobalInputMap
    for mode in 1:3
        for local_input in 1:3
            gid = SY.TimedHybridAutomata.get_global_input_id(gim, mode, local_input)
            @test gid > 0
            typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
            @test typ == :continuous
            @test info[1] == mode
        end
    end
    for tr in 1:6
        gid = SY.TimedHybridAutomata.get_switching_global_id(gim, tr)
        @test gid > 0
        typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
        @test typ == :switching
        @test info == tr
    end

    # 3. Vérification des modèles symboliques dynamiques et temporels
    mode_n_input = [25, 7*7, 25]  
    for (idxsm, sm) in enumerate(hybrid_model.symmodels)
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(sm) == mode_n_input[idxsm]
    end
    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end

    # 4. Vérification de la correspondance états augmentés <-> indices
    for i in 1:length(hybrid_model.int2aug_state)
        @test haskey(hybrid_model.aug_state2int, hybrid_model.int2aug_state[i])
    end

    # 6. Vérification des utilitaires find_symbolic_state, find_time_index, extract_spatial_part, extract_temporal_part
    for mode in 1:3
        sm = hybrid_model.symmodels[mode]
        tm = hybrid_model.time_symbolic_models[mode]
        # Prendre un état concret dans le domaine
        x0 = [0.0, 0.0]
        idx = SY.TimedHybridAutomata.find_symbolic_state(sm, x0)
        @test idx > 0
        t0 = tm.tsteps[1]
        tidx = SY.TimedHybridAutomata.find_time_index(tm, t0)
        @test tidx == 1
        # Cas hors domaine
        idx_bad = SY.TimedHybridAutomata.find_symbolic_state(sm, [100.0, 100.0])
        @test idx_bad == 0
        # extract_spatial_part et extract_temporal_part sur une garde
        guard = guards[mode]
        sp = SY.TimedHybridAutomata.extract_spatial_part(guard)
        tp = SY.TimedHybridAutomata.extract_temporal_part(guard)
        @test sp.lb == guard.lb[1:end-1]
        @test tp == [guard.lb[end], guard.ub[end]]
        # get_time_indices_from_interval
        indices = SY.TimedHybridAutomata.get_time_indices_from_interval(tm, tp)
        @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    end

    

    #  transitions invalides, états hors domaine, etc.
    for mode in 1:3
        sm = hybrid_model.symmodels[mode]
        # find_symbolic_state sur un point très loin
        @test SY.TimedHybridAutomata.find_symbolic_state(sm, [999.0, -999.0]) == 0
    end

    # 9. Vérification exhaustive des entrées globales
    for gid in 1:gim.total_inputs
        typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
        if typ == :continuous
            @test 1 <= info[1] <= 27
            @test 1 <= info[2] <= 99
        elseif typ == :switching
            @test 1 <= info <= 6
        else
            @test typ == :invalid
        end
    end

    # 10. Vérification de la cohérence des indices d'états et de temps dans l'automate
    for i in 1:length(hybrid_model.int2aug_state)
        aug = hybrid_model.int2aug_state[i]
        mode = aug[3]
        t_idx = aug[2]
        @test 1 <= mode <= 3
        @test 1 <= t_idx <= length(hybrid_model.time_symbolic_models[mode].tsteps)
    end
    
end
end