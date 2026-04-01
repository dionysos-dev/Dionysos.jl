export MPPIConfig, MPPIGenerator, get_seed
# je dois (et globalement dans absolument tous mes nouveaux codes remplacer les any[])
# Mettre à jour u_nom = u_new à chaque itération, tout en gardant best_cand séparémen
import Dionysos
import Random

const DI = Dionysos
const ST = DI.System

# ------------------------------------------------------------------
#  MPPI heuristic generator
#
#  Cette V1 implemente volontairement une version tres simple de MPPI.
#
#  Le point cle est le suivant :
#  MPPI est ici un "wrapper" autour d'un autre generateur.
#  Ce generateur seed fournit une premiere sequence de controles via
#  `u_traj`, qui sert de warm-start. En revanche, on ne reutilise pas
#  `x_traj` comme reference physique : toutes les trajectoires MPPI
#  sont reconstruites par simulation concrete avec la dynamique fournie
#  dans la configuration.
# ------------------------------------------------------------------

# Petit helper d'identite pour les cas sans periodicite. On garde cette
# fonction separee pour pouvoir proposer un constructeur pratique de
# `MPPIConfig` sans forcer l'utilisateur a toujours fournir un callback
# `wrap_state`.
_identity_wrap_state(problem, x) = x # je devrais avoir une manière générique de traiter ça

# ------------------------------------------------------------------
#  Configuration MPPI
#
#  Cette configuration reste volontairement compacte.
#  Toute la personnalisation passe par quelques callbacks simples :
#  - comment obtenir l'etat initial concret ;
#  - comment simuler la dynamique discrete ;
#  - comment echantillonner du bruit de controle ;
#  - comment projeter un controle dans le domaine admissible ;
#  - comment evaluer le cout total d'une trajectoire ;
#  - comment decider si la trajectoire est un succes ;
#  - comment enrouler l'etat si le probleme contient des dimensions
#    periodiques.
#
#  Ce choix garde la V1 generique sans multiplier les couches
#  d'abstraction inutiles.
# ------------------------------------------------------------------
struct MPPIConfig{FX0, FDYN, FNOISE, FPROJ, FCOST, FSUCC, FWRAP}
    # Pas de temps associe a la trajectoire candidate renvoyee.
    Δt::Float64

    # Horizon d'optimisation exprime en nombre de controles.
    nstep::Int # -> jouer avec l'horizon à travers la cost function

    # Nombre de rollouts stochastiques par iteration MPPI.
    nsamples::Int # boucle 1

    # Nombre d'iterations MPPI.
    niter::Int # boucle 2

    # Temperature MPPI.
    # Plus λ est petit, plus l'algorithme favorise agressivement les samples de faible cout.
    λ::Float64

    # x0_provider(problem) -> x0
    # Permet de definir l'etat initial 
    x0_provider::FX0

    # discrete_dynamics(problem, x, u, k, Δt) -> x_next (ça pourrait etre différent de celle obtenu pour la traj nominal, par exemple plus riche)
    discrete_dynamics::FDYN

    # noise_sampler(rng, u, k) -> ϵ ; ϵ ~ N(0,Σ) (dans les faits on utilisera pas d'autre distribution)
    # Fournit une perturbation de controle pour le pas k.
    noise_sampler::FNOISE

    # project_input(u) -> u_proj (forcer le control à rester admissible)
    # Utilise pour re-clamper / re-projeter les controles apres
    # perturbation et apres mise a jour MPPI.
    project_input::FPROJ

    # trajectory_cost(problem, cand) -> J
    # on laisse l'utilisateur definir librement un cout total de trajectoire. 
    trajectory_cost::FCOST

    # success_fun(problem, cand) -> Bool
    success_fun::FSUCC

    # wrap_state(problem, x) -> x_wrapped
    # Ce callback permet de normaliser l'etat apres chaque pas de
    # simulation, par exemple pour les dimensions periodiques.
    wrap_state::FWRAP
end

function MPPIConfig(
    Δt::Real,
    nstep::Integer,
    nsamples::Integer,
    niter::Integer,
    λ::Real,
    x0_provider,
    discrete_dynamics,
    noise_sampler,
    project_input,
    trajectory_cost,
    success_fun,
    wrap_state = _identity_wrap_state,
)
    return MPPIConfig{
        typeof(x0_provider), # je dois changer les immondes typeof
        typeof(discrete_dynamics),
        typeof(noise_sampler),
        typeof(project_input),
        typeof(trajectory_cost),
        typeof(success_fun),
        typeof(wrap_state),
    }(
        Float64(Δt),
        Int(nstep),
        Int(nsamples),
        Int(niter),
        Float64(λ),
        x0_provider,
        discrete_dynamics,
        noise_sampler,
        project_input,
        trajectory_cost,
        success_fun,
        wrap_state,
    )
end

# ------------------------------------------------------------------
#  Generateur MPPI
#
#  Le generateur stocke :
#  - le probleme courant ;
#  - la configuration MPPI ;
#  - un generateur seed externe ;
#  - la candidate finale ;
#  - un booleen de succes ;
#  - le temps de calcul ;
#  - un petit bloc de diagnostics.
# ------------------------------------------------------------------
mutable struct MPPIGenerator{P, C, SG, D} <: AbstractHeuristicGenerator
    # Probleme concret a resoudre.
    problem::Union{Nothing, P}

    # Configuration de l'algorithme MPPI.
    config::C

    # Generateur seed utilise pour produire la sequence de controles
    # initiale servant de warm-start.
    seed_generator::SG

    # Trajectoire brute renvoyee par le generateur seed.
    seed::Union{Nothing, CandidateTrajectory}

    # Candidate finale renvoyee au reste du pipeline.
    candidate::Union{Nothing, CandidateTrajectory}

    # Vrai si la trajectoire finale satisfait `cfg.success_fun`.
    success::Bool

    # Temps total passe dans `generate!`.
    solve_time_sec::Float64

    # Petit conteneur de diagnostics.
    # On le garde tres simple pour cette V1.
    diagnostics::D
end

function MPPIGenerator(problem, config, seed_generator; diagnostics = (;))
    return MPPIGenerator{Any, typeof(config), typeof(seed_generator), Any}(
        problem,
        config,
        seed_generator,
        nothing,
        nothing,
        false,
        0.0,
        diagnostics,
    )
end

# ------------------------------------------------------------------
#  Ajustement de l'horizon
#
#  Cette V1 adopte une politique volontairement simple :
#  - si la seed est trop longue, on la tronque ;
#  - si elle est trop courte, on repete le dernier controle ;
#  - sinon, on la garde telle quelle.
#
#  On n'introduit ni interpolation ni heuristique plus fine, afin de
#  garder un comportement lisible et previsible.
# ------------------------------------------------------------------
function _pad_or_trim_controls(u_seq, nstep::Int) # je devrais peut etre avoir exactement la meme taille que ma candidate traj mais enfaite je préfère gérer ça avec la cost function, ça donne un degré de liberté en plus
    length(u_seq) == nstep && return collect(u_seq)
    length(u_seq) > nstep && return collect(u_seq[1:nstep])

    padded = collect(u_seq)
    last_u = padded[end]
    while length(padded) < nstep
        push!(padded, last_u)
    end
    return padded
end

# ------------------------------------------------------------------
#  Rollout concret d'une sequence de controles
#
#  Ce choix est important :
#  meme si la seed provient d'une trajectoire abstraite, les trajectoires
#  evaluees par MPPI sont toujours des trajectoires concretes obtenues
#  par simulation.
# ------------------------------------------------------------------
function _rollout_controls(problem, cfg::MPPIConfig, x0, u_seq)
    x = cfg.wrap_state(problem, x0)

    xs = Vector{typeof(x)}(undef, 0)
    push!(xs, x)

    us = Any[]

    for (k, u) in enumerate(u_seq)
        push!(us, u)
        x = cfg.discrete_dynamics(problem, x, u, k, cfg.Δt)
        x = cfg.wrap_state(problem, x)
        push!(xs, x)
    end

    return CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = cfg.Δt,
        source = :mppi,
        metadata = (
            ;
            nstep = cfg.nstep,
            nsamples = cfg.nsamples,
            niter = cfg.niter,
        ),
    )
end

function _truncate_candidate_at_first_target_hit(problem, cfg::MPPIConfig, cand)
    cand === nothing && return nothing
    hasproperty(problem, :target_set) || return cand

    xs = collect(ST.enum_elems(cand.x_traj))
    hit_idx = findfirst(x -> (cfg.wrap_state(problem, x) ∈ problem.target_set), xs)

    if hit_idx === nothing || hit_idx <= 1 || hit_idx >= length(xs)
        return cand
    end

    us = collect(ST.enum_elems(cand.u_traj))

    return CandidateTrajectory(
        ST.Trajectory(xs[1:hit_idx]),
        ST.Trajectory(us[1:(hit_idx - 1)]);
        Ts = cand.Ts,
        source = cand.source,
        metadata = cand.metadata,
    )
end

# ------------------------------------------------------------------
#  Une iteration MPPI
#
#  Cette fonction implemente une mise a jour tres simple :
#  1. on part de la sequence nominale `u_nom` ;
#  2. on echantillonne `nsamples` sequences perturbees ;
#  3. on simule chaque trajectoire concrete ;
#  4. on evalue son cout total ;
#  5. on construit les poids exponentiels classiques de MPPI ;
#  6. on fait une moyenne ponderee des perturbations ;
#  7. on applique cette correction a la sequence nominale.
#
# ------------------------------------------------------------------
function _mppi_update(problem, cfg::MPPIConfig, x0, u_nom, rng)
    nsamples = cfg.nsamples
    horizon = length(u_nom)

    eps_samples = Vector{Vector{Any}}(undef, nsamples)
    costs = Vector{Float64}(undef, nsamples)

    for s in 1:nsamples # on peut multithread cette boucle (et on va)
        eps_seq = Any[]
        u_roll = Any[]

        for k in 1:horizon
            ϵ = cfg.noise_sampler(rng, u_nom[k], k)
            push!(eps_seq, ϵ)
            push!(u_roll, cfg.project_input(u_nom[k] + ϵ))
        end

        cand = _rollout_controls(problem, cfg, x0, u_roll)
        costs[s] = Float64(cfg.trajectory_cost(problem, cand))
        eps_samples[s] = eps_seq
    end

    β = minimum(costs) # on devrait peut etre ajouter une safeguard numérique
    weights = exp.(-(costs .- β) ./ cfg.λ)
    weights ./= sum(weights)

    u_new = Any[]

    for k in 1:horizon
        # en multipliant la premiere perturbation par 0.0.
        # Cela reste simple et fonctionne bien pour les types usuels
        # attendus ici : scalaires, vecteurs statiques, etc.
        δu = 0.0 * eps_samples[1][k]

        for s in 1:nsamples
            δu += weights[s] * eps_samples[s][k]
        end

        push!(u_new, cfg.project_input(u_nom[k] + δu))
    end

    return u_new, (; sample_best_cost = β)
end

function set_problem!(gen::MPPIGenerator, problem)
    gen.problem = problem
    set_problem!(gen.seed_generator, problem)
    gen.seed = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    gen.diagnostics = (;)
    return gen
end

# ------------------------------------------------------------------
#  1. reinitialiser l'etat du generateur ;
#  2. lancer le generateur seed ;
#  3. recuperer sa sequence de controles comme warm-start ; # à terme on devrait extend : on devrait pouvoir récupérer pluseorus sequences de control et faire du multiple MPPI (warm-start), mais vu les résultats actuelle ça me semble pas si fou
#  4. reconstruire une premiere trajectoire concrete ;
#  5. iterer MPPI ;
#  6. conserver la meilleure candidate ;
#  7. definir le succes uniquement a partir de la candidate finale.
#
#  Point important :
#  on n'exige pas que le seed_generator ait `success == true`.
#  Tant qu'il fournit une trajectoire seed non nulle, son `u_traj`
#  peut etre utile comme point de depart pour MPPI.
# ------------------------------------------------------------------
function generate!(gen::MPPIGenerator)
    cfg = gen.config

    # def
    @assert gen.problem !== nothing "Call set_problem!(gen, problem) first."
    @assert cfg.Δt > 0.0
    @assert cfg.nstep >= 1
    @assert cfg.nsamples >= 1
    @assert cfg.niter >= 1
    @assert cfg.λ > 0.0

    gen.seed = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    gen.diagnostics = (;)

    t0 = time()

    # On commence par produire la seed. Cette seed sert uniquement de
    # warm-start en controle. Elle n'a pas besoin d'etre deja "reussie".
    generate!(gen.seed_generator)
    seed = get_trajectory(gen.seed_generator)
    gen.seed = seed

    # Sans seed, cette V1 ne lance pas MPPI. Le contrat reste simple :
    # pas de trajectoire seed -> pas de sequence nominale initiale. # à terme on pourrait mettre ça en option, dans les articles il précise que l'on peut tout à fait init MPPI avec une u_seq constante
    if seed === nothing
        gen.solve_time_sec = time() - t0
        gen.diagnostics = (; seed_available = false)
        return gen
    end

    # On recupere l'etat initial concret et la sequence nominale issue du
    # seed, puis on l'ajuste a l'horizon demande.
    problem = gen.problem
    x0 = cfg.x0_provider(problem)
    u_nom = collect(ST.enum_elems(seed.u_traj))
    u_nom = _pad_or_trim_controls(u_nom, cfg.nstep)

    # Premiere candidate concrete : on re-simule toujours la trajectoire
    # a partir de la dynamique concrete, meme si le seed venait d'un mode abstrait.
    best_cand = _rollout_controls(problem, cfg, x0, u_nom)
    seed_cost = Float64(cfg.trajectory_cost(problem, best_cand))
    best_cost = seed_cost

    rng = Random.default_rng()
    iterations_done = 0

    for it in 1:cfg.niter
        iterations_done = it

        # Une iteration MPPI produit une nouvelle sequence nominale
        # candidate, calculee a partir d'une moyenne ponderee des
        # perturbations echantillonnees.
        u_new, _ = _mppi_update(problem, cfg, x0, u_nom, rng)

        cand = _rollout_controls(problem, cfg, x0, u_new)
        cost = Float64(cfg.trajectory_cost(problem, cand))

        # version défensive (1) 

        # Pour cette V1, on ne remplace le nominal courant que si la
        # nouvelle candidate ameliore effectivement le meilleur cout vu
        # jusqu'ici.
        if cost < best_cost # théoriquement ça n'a pas trop de sens, je devrais tester si c'est performant ou pas
            best_cost = cost
            best_cand = cand
            u_nom = u_new
        end

        # version classique MPPI (2)
        # u_nom = u_new
        # if cost < best_cost 
        #     best_cost = cost
        #     best_cand = cand
        # end

        # Si la meilleure candidate satisfait deja le critere de succes,
        # on s'arrete. Cela suffit pour une premiere implementation.
        if cfg.success_fun(problem, best_cand) # théoriquement ça n'a pas trop de sens aussi
            break
        end
    end

    final_cand = _truncate_candidate_at_first_target_hit(problem, cfg, best_cand)

    gen.candidate = final_cand
    gen.success = cfg.success_fun(problem, final_cand)
    gen.solve_time_sec = time() - t0
    gen.diagnostics = (
        ;
        seed_cost = seed_cost,
        best_cost = best_cost,
        final_cost = best_cost,
        iterations = iterations_done,
    )
    return gen
end

get_trajectory(gen::MPPIGenerator) = gen.candidate
get_seed(gen::MPPIGenerator) = gen.seed
get_success(gen::MPPIGenerator) = gen.success
get_solve_time(gen::MPPIGenerator) = gen.solve_time_sec
