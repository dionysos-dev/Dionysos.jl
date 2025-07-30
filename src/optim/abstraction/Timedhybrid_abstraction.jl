export TemporalHybridSymbolicModelAbstraction

module TemporalHybridSymbolicModelAbstraction

import HybridSystems
import HybridSystems: HybridSystem
import MathematicalSystems
import StaticArrays: SVector
import Dionysos

struct ProblemSpecs{F}
    initial_state::Tuple{AbstractVector{Float64}, Float64, Int} # état augmenté ([x], t, mode_id)
    Xs_target::Vector{<:Dionysos.Utils.HyperRectangle} # must be changed to admit more complex target sets
    Ts_target::Vector{<:Dionysos.Utils.HyperRectangle}
    Ns_target::Vector{Int}
    concret_cost_fun::F
end

# =================================
# SYNTHÈSE DE CONTRÔLEUR
# =================================

"""Construit le problème de contrôle concret"""
function build_concrete_problem(Concret_specifications::ProblemSpecs,) 
    concrete_initial_set = Concret_specifications.initial_state
    concrete_target_set = (Concret_specifications.Xs_target, Concret_specifications.Ts_target, Concret_specifications.Ns_target)
    concret_cost_fun = Concret_specifications.concret_cost_fun

    return Dionysos.Problem.OptimalControlProblem(
        Concret_specifications,
        concrete_initial_set,
        concrete_target_set,
        nothing,
        concret_cost_fun,
        Dionysos.Problem.Infinity(),
    )
end

"""Construit le problème de contrôle abstrait"""
function build_abstract_problem(concrete_problem::Dionysos.Problem.OptimalControlProblem, symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel) 
    aug_state = concrete_problem.initial_set
    abstract_initial_set = [Dionysos.Symbolic.TimedHybridAutomata.get_abstract_state(symmodel, aug_state)] 
    target_set = concrete_problem.target_set 
    abstract_target_set = Dionysos.Symbolic.TimedHybridAutomata.get_states_from_set(symmodel, target_set[1], target_set[2], target_set[3]) 

    return Dionysos.Problem.OptimalControlProblem(
        symmodel,
        abstract_initial_set,
        abstract_target_set,
        concrete_problem.state_cost,
        build_abstract_transition_cost(symmodel, concrete_problem.transition_cost), 
        concrete_problem.time,
    ) 
end

"""Résout le problème abstrait et calcule l'ensemble contrôlable"""
function solve_abstract_problem(abstract_problem)

    abstract_controller, controllable_set_symbols, _, value_per_node = Dionysos.Symbolic.compute_worst_case_cost_controller(
        abstract_problem.system.autom,
        abstract_problem.target_set;
        initial_set = abstract_problem.initial_set,
        sparse_input = false, 
        cost_function = abstract_problem.transition_cost,
    )
    
    println("value init_state : ",value_per_node[abstract_problem.initial_set[1]])

    if ⊆(abstract_problem.initial_set, controllable_set_symbols)
        println("✅ Problem is solvable: initial set is controllable")
    else
        println("⚠️ Warning: initial set is only partially controllable")
    end
    
    return abstract_controller, controllable_set_symbols
end

"""Génère le contrôleur concret à partir du contrôleur abstrait"""
function solve_concrete_problem(symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel, abstract_controller)
    
    function concrete_controller(aug_state)
        (x, t, k) = aug_state
        abstract_aug_state = Dionysos.Symbolic.TimedHybridAutomata.get_abstract_state(symmodel, aug_state)
        
        if Dionysos.System.is_defined(abstract_controller, abstract_aug_state) == false
            println("⚠️ No action available for state $aug_state")
            return nothing
        end

        abstract_input = Dionysos.System.get_control(abstract_controller, abstract_aug_state)

        if Dionysos.Symbolic.TimedHybridAutomata.is_switching_input(symmodel.global_input_map, abstract_input)
            transition_id = symmodel.global_input_map.global_to_switching[abstract_input]
            label = symmodel.global_input_map.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            return label
        else
            return Dionysos.Symbolic.TimedHybridAutomata.get_concrete_input(symmodel, abstract_input, k)
        end
    end
    
    isdef = (aug_state) -> true # for the moment we assume the controller is defined for all augmented states
    return Dionysos.System.BlackBoxContinuousController(concrete_controller, isdef)
end

# =================================
# FONCTIONS DE COÛT ABSTRAIT
# =================================

function build_abstract_transition_cost(symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel, concrete_cost_fun)
    function abstract_transition_cost(state_int, input_int)
        (x,t,k) = Dionysos.Symbolic.TimedHybridAutomata.get_concrete_state(symmodel, state_int)
        aug_state = (x, t, k) 
        if Dionysos.Symbolic.TimedHybridAutomata.is_switching_input(symmodel.global_input_map, input_int)
            transition_id = symmodel.global_input_map.global_to_switching[input_int]
            label = symmodel.global_input_map.switch_labels[transition_id] # eg. "SWITCH source_mode_id -> target_mode_id"
            u = label
        else
            u = Dionysos.Symbolic.TimedHybridAutomata.get_concrete_input(symmodel, input_int, k)
        end
        return concrete_cost_fun(aug_state, u)
    end
    return abstract_transition_cost
end

# =================================
# FONCTION PRINCIPALE DE RÉSOLUTION
# =================================

"""
Résout le problème de contrôle hybride temporel complet.

# Arguments
- `tasks` : Liste des tâches
- `symmodel` : Modèle symbolique de base
- `tstep` : Pas de discrétisation temporelle

# Retourne
Contrôleur concret pour états augmentés
"""
function solve(hs::HybridSystem, growth_bounds::SVector{}, param_discretization, concret_specs::ProblemSpecs) 
    #println("🔄 Construction du modèle temporel...")
    hybrid_symmodel = Dionysos.Symbolic.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretization,
    )

    # Problème concret
    #println("🔄 Construction du problème de contrôle concret...")
    concrete_problem = build_concrete_problem(concret_specs)  
    
    # Problème abstrait
    #println("🔄 Construction du problème de contrôle abstrait...")
    abstract_problem = build_abstract_problem(concrete_problem, hybrid_symmodel) 
    
    
    # Résolution abstraite
    #println("🔄 Résolution du problème de contrôle abstrait...")
    abstract_controller, _ = solve_abstract_problem(abstract_problem) 
    
    # Contrôleur concret
    #println("🔄 Génération du contrôleur concret...")
    concrete_controller = solve_concrete_problem(hybrid_symmodel, abstract_controller) 
    
    return concrete_controller
end

# =================================
# SIMULATION EN BOUCLE FERMÉE
# =================================

"""Calcule l'état augmenté suivant"""
function get_next_aug_state(hs::HybridSystem, aug_state, u, tm, map_sys)
    (x, t, k) = aug_state
    #println("\n \nCurrent augmented state: ", aug_state)
    #println("Current input: ", u)

    if isa(u, AbstractString) && occursin("SWITCH", u)
        #println("Switching input detected: ", u)
        # Extraire l'id de transition à partir du label
        m = match(r"SWITCH (\d+) -> (\d+)", String(u))
        if m === nothing
            error("Format de label SWITCH non reconnu: $u")
        end
        source_mode = parse(Int, m.captures[1])
        target_mode = parse(Int, m.captures[2])

        # Trouver l'index de la transition correspondante
        transitions = collect(HybridSystems.transitions(hs.automaton))
        transition_id = findfirst(tr -> HybridSystems.source(hs.automaton, tr) == source_mode &&
                                         HybridSystems.target(hs.automaton, tr) == target_mode, transitions)
        if transition_id === nothing
            error("Transition non trouvée pour $u")
        end
        transition = transitions[transition_id]

        # Appliquer la reset map de la transition
        reset_map = HybridSystems.resetmap(hs, transition)
        augmented_source_state = vcat(x, t)
        reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
        next_x = reset_result[1:end-1]
        next_t = reset_result[end]
        next_k = target_mode
        return (next_x, next_t, next_k)
    else
        # Évolution continue
        next_t = tm.is_active ? t + tm.tstep : 0.0
        next_x = map_sys(x, u, tm.tstep)
        return (next_x, next_t, k)
    end
end

"""
Génère une trajectoire en boucle fermée.

# Arguments
- `tasks` : Liste des tâches
- `continuous_time_system` : Système dynamique continu
- `controller` : Contrôleur synthétisé
- `aug_state_0` : État initial augmenté
- `nstep` : Nombre maximal d'étapes
- `stopping` : Fonction d'arrêt

# Retourne
Tuple (trajectoire_états, trajectoire_commandes)
"""
function get_closed_loop_trajectory(
    symmodel::Dionysos.Symbolic.TimedHybridAutomata.TemporalHybridSymbolicModel,
    hs::HybridSystem,
    problem_specs::ProblemSpecs,
    controller,
    aug_state_0,
    nstep;
    stopping = (x) -> false,
)
    #println("\nwe are in get_closed_loop_trajectory\n")
    aug_state_traj, u_traj = [aug_state_0], []
    aug_state = aug_state_0
    
    nmodes = HybridSystems.nmodes(hs.automaton)
    dynamics = [HybridSystems.mode(hs, k).systems[1].f for k in 1:nmodes]
    maps_sys = [Dionysos.System.simulate_control_map(dynamics[i]) for i in 1:nmodes]


    for _ in 1:nstep
        #println("stopping works ? ", stopping(problem_specs, aug_state))
        stopping(problem_specs, aug_state) && break
        

        #println("Current augmented state: ", aug_state)
        u = controller.f(aug_state) 
        #println("control works ? ", controller.f(aug_state) )

        u === nothing && break
        (_,_,k) = aug_state
        aug_state = get_next_aug_state(hs, aug_state, u, symmodel.time_symbolic_models[k], maps_sys[k])
      
        push!(aug_state_traj, aug_state)
        push!(u_traj, u)
    end
    
    return aug_state_traj, u_traj
end

"""Vérifie si l'objectif final est atteint
"""
function reached(specs::ProblemSpecs, aug_state)
    (x, t, k) = aug_state
    # Vérifie si le mode courant fait partie des modes cibles
    idx = findfirst(==(k), specs.Ns_target)
    if isnothing(idx)
        return false
    end
    # Vérifie inclusion spatiale et temporelle
    X_target = specs.Xs_target[idx]
    T_target = specs.Ts_target[idx]
    in_X = x ∈ X_target
    in_T = t ≥ T_target.lb[1] && t ≤ T_target.ub[1]
    return in_X && in_T
end

end
