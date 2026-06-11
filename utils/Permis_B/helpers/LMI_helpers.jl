# ----------------------------------------------------------------------
# Helpers LMI (mode symbolique uniquement)
# ----------------------------------------------------------------------

#=
Guide detaille - Helpers LMI symboliques
========================================

Ce fichier implemente uniquement la voie LMI symbolique.

Contexte d'utilisation
----------------------
Ce fichier est charge par `permis_B.jl` et suppose les alias deja definis:
- `UT` (Dionysos.Utils)
- `ST` (Dionysos.System)
- `SY` (Dionysos.Symbolic)

Probleme traite
---------------
On dispose d'une trajectoire nominale `(x[k], u[k])` et d'une ellipsoide terminale `E_N`.
On cherche ensuite des ellipsoides `E_k` en remontant le temps, telles que chaque transition
`E_k --u[k]--> E_{k+1}` soit certifiee par une LMI.

Etapes d'une transition k
-------------------------
1) Lineariser localement le modele symbolique autour de `(xk, uk)` sur un voisinage
   intervalle (`ΔX`, `ΔU`, `ΔW`).
2) Resoudre `SY.transition_backward(...)` pour trouver `E_prev` et un controle local `kappa`.
3) Si succes: remplacer `E_next` par `E_prev` et continuer. Sinon: stopper la chaine.

Parametres critiques
--------------------
- `ΔX`, `ΔU`, `ΔW`: taille des voisinages de linearisation.
- `maxδx`, `maxδu`: bornes internes utilisees dans la LMI.
- `λ`: regularisation numerique du SDP.

Causes frequentes d'infeasibilite
---------------------------------
- voisinages trop larges (`ΔX`, `ΔU`),
- regularisation insuffisante (`λ` trop faible),
- incoherence entre trajectoire nominale et modele symbolique local.
=#

using JuMP
using Clarabel
import IntervalArithmetic as IA
import MathOptInterface as MOI
import LinearAlgebra as LA
import MathematicalSystems as MS
import Random

"""
    LMIContext

Contexte complet pour la synthese backward:
- dynamique symbolique (`fsymbolic`, `x`, `u`, `w`),
- voisinages de linearisation (`ΔX`, `ΔU`, `ΔW`),
- formats des contraintes (`Uformat`, `Wformat`),
- cout quadratique (`S`), solveur SDP (`sdp_opt`),
- hyperparametres (`maxδx`, `maxδu`, `λ`).
"""
Base.@kwdef struct LMIContext
    fsymbolic::Any
    x::Any
    u::Any
    w::Any
    ΔX::Any
    ΔU::Any
    ΔW::Any
    Uformat::Any
    Wformat::Any
    S::Any
    sdp_opt::Any
    maxδx::Float64 = 100.0
    maxδu::Float64 = 100.0
    λ::Float64 = 0.01
    use_noise::Bool = true
    use_logdet::Bool = true
    trace_reward::Float64 = 0.0
end

"""
    construire_contexte_lmi(; ...)

Cree un `LMIContext` pret a l'emploi.
`transition_cost` accepte:
- soit une matrice deja construite,
- soit une fonction quadratique Dionysos (convertie via `UT.get_full_psd_matrix`).
"""
function construire_contexte_lmi(;
    fsymbolic,
    x,
    u,
    w,
    ΔX,
    ΔU,
    ΔW,
    Uformat,
    Wformat,
    transition_cost,
    sdp_opt = nothing,
    maxδx = 100.0,
    maxδu = 20.0,
    λ = 0.01,
)
    S =
        transition_cost isa AbstractMatrix ? transition_cost :
        UT.get_full_psd_matrix(transition_cost)

    return LMIContext(;
        fsymbolic = fsymbolic,
        x = x,
        u = u,
        w = w,
        ΔX = ΔX,
        ΔU = ΔU,
        ΔW = ΔW,
        Uformat = Uformat,
        Wformat = Wformat,
        S = S,
        sdp_opt = sdp_opt,
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
    )
end

"""
    matrice_cout_identite(nx, nu)

Construit la matrice de cout:
`||x||^2 + ||u||^2 + 1`.
"""
function matrice_cout_identite(nx::Int, nu::Int)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

"""
    lineariser_localement(ctx, xnew; unew, wnew)

Linearise le modele symbolique autour de `(xnew, unew, wnew)`:
- construit `Xbar/Ubar/Wbar`,
- appelle `ST.buildAffineApproximation`,
- retourne `(affineSys, L, uvec, wvec)`.
"""
function lineariser_localement(
    ctx::LMIContext,
    xnew;
    unew = zeros(length(ctx.u)),
    wnew = zeros(length(ctx.w)),
)
    xvec = collect(xnew)
    uvec = collect(unew)
    wvec = collect(wnew)

    Xbar = IA.IntervalBox(xvec .+ ctx.ΔX)
    Ubar = IA.IntervalBox(uvec .+ ctx.ΔU)
    Wbar = IA.IntervalBox(wvec .+ ctx.ΔW)

    affineSys, L = ST.buildAffineApproximation(
        ctx.fsymbolic,
        ctx.x,
        ctx.u,
        ctx.w,
        xvec,
        uvec,
        wvec,
        Xbar,
        Ubar,
        Wbar,
        #remainder_model = remainder_model,
    )

    return affineSys, L, uvec, wvec
end

"""
    _eval_local_controller(kappa, x)

Evalue un controleur local sur l'etat `x`.
Prend en charge:
- `MS.AffineMap` (champs `A` + `b` ou `A` + `c`),
- matrice augmentee `[K  ell]`,
- fallback via `ST.get_control`.
"""
function _eval_local_controller(kappa, x::AbstractVector)
    xvec = collect(Float64, x)

    if hasproperty(kappa, :A) && hasproperty(kappa, :b)
        A = Matrix{Float64}(getproperty(kappa, :A))
        b = vec(Float64.(getproperty(kappa, :b)))
        return collect(A * xvec + b)
    end
    if hasproperty(kappa, :A) && hasproperty(kappa, :c)
        A = Matrix{Float64}(getproperty(kappa, :A))
        b = vec(Float64.(getproperty(kappa, :c)))
        return collect(A * xvec + b)
    end
    if kappa isa AbstractMatrix
        K = Matrix{Float64}(kappa[:, 1:(end - 1)])
        ell = vec(Float64.(kappa[:, end]))
        return collect(K * xvec + ell)
    end

    return collect(Float64.(ST.get_control(kappa, xvec)))
end

"""
    creer_ellipsoide_terminal(state_list; rayon, P)

Cree l'ellipsoide terminale centree sur `state_list[end]`.
Par defaut: boule isotrope de rayon `rayon`.
"""
function creer_ellipsoide_terminal(state_list; rayon = 0.25, P = nothing)
    isempty(state_list) && error("state_list est vide.")

    c = collect(state_list[end])
    nx = length(c)
    Pmat =
        P === nothing ? Matrix{Float64}(LA.I, nx, nx) * (1.0 / rayon^2) : Matrix{Float64}(P)

    return UT.Ellipsoid(Pmat, c)
end

"""
    synthetiser_transitions_backward(state_list, input_list, ctx; ...)

Remonte la trajectoire nominale en validant les transitions backward par LMI.
Renvoie un NamedTuple avec succes/echec, indices, couts et ellipsoides.
"""
function synthetiser_transitions_backward(
    state_list,
    input_list,
    ctx::LMIContext;
    E_terminal = nothing,
    stride::Int = 1,
    verbose::Bool = true,
)
    n = length(state_list)
    n >= 2 || error("Il faut au moins deux etats dans state_list.")
    stride >= 1 || error("stride doit etre >= 1.")

    E_next = E_terminal === nothing ? creer_ellipsoide_terminal(state_list) : E_terminal

    ellipsoides = UT.Ellipsoid[E_next]
    controlleurs = Any[]
    couts = Float64[]
    indices = Int[]

    for k in (n - 1):(-stride):1
        xk = collect(state_list[k])
        uk = k <= length(input_list) ? collect(input_list[k]) : zeros(length(ctx.u))

        affineSys, L, uvec, _ = lineariser_localement(ctx, xk; unew = uk)
        verbose && println("k=", k, " | Lips=", L)

        E_try, kappa_try, cout_try = SY.transition_backward(
            affineSys,
            E_next,
            xk,
            uvec,
            ctx.Uformat,
            ctx.Wformat,
            ctx.S,
            L,
            ctx.sdp_opt;
            maxδx = ctx.maxδx,
            maxδu = ctx.maxδu,
            λ = ctx.λ,
        )

        if E_try === nothing || kappa_try === nothing
            verbose && println("Transition impossible pour k=", k, " (x[k+1] -> x[k]).")
            return (;
                success = false,
                failed_k = k,
                ellipsoides = ellipsoides,
                controlleurs = controlleurs,
                couts = couts,
                indices = indices,
            )
        end

        push!(ellipsoides, E_try)
        push!(controlleurs, kappa_try)
        push!(couts, cout_try)
        push!(indices, k)
        E_next = E_try

        println("the volume is ", UT.get_volume(E_next), "\n\n")

        verbose && println("Transition OK k=", k, ", cout=", cout_try)
    end

    return (;
        success = true,
        failed_k = nothing,
        ellipsoides = ellipsoides,
        controlleurs = controlleurs,
        couts = couts,
        indices = indices,
    )
end

"""
    _sample_points_uniform_in_ellipsoid(E, n; rng = Random.default_rng())

Echantillonne `n` points (quasi-uniformement) dans l'ellipsoide
`E = {x : (x-c)'P(x-c) <= 1}`.
"""
function _sample_points_uniform_in_ellipsoid(
    E::UT.Ellipsoid,
    n::Int;
    rng = Random.default_rng(),
)
    n >= 1 || error("n doit etre >= 1.")
    nx = length(E.c)
    L = LA.cholesky(LA.Symmetric(Matrix{Float64}(E.P))).L
    c = vec(Float64.(E.c))

    pts = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        v = Random.randn(rng, nx)
        nv = LA.norm(v)
        while nv <= 1e-14
            v = Random.randn(rng, nx)
            nv = LA.norm(v)
        end
        # Rayon pour une loi uniforme dans la boule unite en dimension nx.
        ρ = Random.rand(rng)^(1 / nx)
        z = (ρ / nv) .* v
        pts[i] = collect(c .+ L \ z)
    end
    return pts
end

"""
    _project_on_input_format(u, Uformat; tol = 1e-9)

Projection radiale de `u` sur l'ensemble admissible défini par `Uformat`,
avec contraintes de type `||U[i] * u||₂ ≤ 1`.

Remarque: ce n'est pas la projection euclidienne exacte, mais une méthode
simple et robuste qui conserve la direction de `u`.
"""
function _project_on_input_format(u::AbstractVector, Uformat; tol::Float64 = 1e-9)
    isempty(Uformat) && return collect(Float64, u)
    uvec = collect(Float64, u)

    max_violation = 0.0
    for Ui in Uformat
        max_violation = max(max_violation, LA.norm(Ui * uvec))
    end

    if max_violation <= 1.0 + tol
        return uvec
    end
    return uvec / max_violation
end

"""
    _build_kappa_and_ellipsoid_maps(transitions)

Construit:
- `k_sorted`: indices temporels k tries croissants,
- `kappa_by_k[k] = κ_k`,
- `E_by_k[k] = E_k` et `E_by_k[k_max+1] = E_terminal`.
"""
function _build_kappa_and_ellipsoid_maps(transitions)
    isempty(transitions.indices) && error("Aucune transition disponible.")
    length(transitions.indices) == length(transitions.controlleurs) ||
        error("Incoherence: indices et controlleurs n'ont pas la meme longueur.")
    length(transitions.ellipsoides) == length(transitions.indices) + 1 ||
        error("Incoherence: il faut |ellipsoides| = |indices| + 1.")

    kappa_by_k = Dict{Int, Any}()
    for (k, κ) in zip(transitions.indices, transitions.controlleurs)
        kappa_by_k[k] = κ
    end
    k_sorted = sort(collect(keys(kappa_by_k)))

    # Convention de stockage de la synthese backward:
    # ellipsoides[1] = E_{kmax+1} (terminale), puis ellipsoides[i+1] = E_{indices[i]}.
    E_by_k = Dict{Int, UT.Ellipsoid}()
    for i in eachindex(transitions.indices)
        E_by_k[transitions.indices[i]] = transitions.ellipsoides[i + 1]
    end
    E_by_k[maximum(k_sorted) + 1] = transitions.ellipsoides[1]

    return k_sorted, kappa_by_k, E_by_k
end

"""
    simuler_kappa_sur_modele_concret(run_result, transitions, ctx; ...)

Verification empirique des lois de controle `κ_k` synthétisées:
1) echantillonne `n_samples` points dans l'ellipsoide initial `E_{k_min}`,
2) applique la chaine de controleurs locaux `κ_k` sur le modele concret discretise,
3) retourne les trajectoires et des statistiques de succes.

Definition d'un succes rollout:
- l'etat reste dans les ellipsoides de la chaine (`inside_chain`),
- l'etat final appartient a l'ellipsoide terminal (`final_in_target`),
- (optionnel) l'etat reste dans le domaine concret (`left_domain == false`).
"""
function simuler_kappa_sur_modele_concret(
    run_result,
    transitions,
    ctx::LMIContext;
    n_samples::Int = 100,
    num_substeps::Int = 5,
    project_inputs::Bool = true,
    check_domain::Bool = true,
    rng = Random.default_rng(),
    wrap_state::Function = identity,
    verbose::Bool = true,
)
    n_samples >= 1 || error("n_samples doit etre >= 1.")
    num_substeps >= 1 || error("num_substeps doit etre >= 1.")

    k_sorted, kappa_by_k, E_by_k = _build_kappa_and_ellipsoid_maps(transitions)
    k_min = minimum(k_sorted)
    k_max = maximum(k_sorted)
    E_init = E_by_k[k_min]
    E_target = E_by_k[k_max + 1]

    x0_samples = _sample_points_uniform_in_ellipsoid(E_init, n_samples; rng = rng)
    disc = ST.discretize_continuous_system(
        run_result.concrete_system,
        run_result.Δt;
        num_substeps = num_substeps,
    )
    fdisc = MS.mapping(disc)

    x_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    u_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    rollout_stats = Vector{NamedTuple}(undef, n_samples)

    for s in 1:n_samples
        x = copy(x0_samples[s])
        x_hist = Vector{Vector{Float64}}([copy(x)])
        u_hist = Vector{Vector{Float64}}()

        inside_chain = true
        left_domain = false

        # Verifie l'appartenance de l'etat initial a E_{k_min}.
        inside_chain &= (x ∈ E_init)

        for k in k_sorted
            if haskey(E_by_k, k)
                inside_chain &= (x ∈ E_by_k[k])
            end

            x_ctrl = wrap_state(x)
            u = _eval_local_controller(kappa_by_k[k], x_ctrl)
            if project_inputs
                u = _project_on_input_format(u, ctx.Uformat)
            end
            push!(u_hist, copy(u))

            x_next = collect(Float64.(fdisc(x, u)))
            push!(x_hist, copy(x_next))

            if haskey(E_by_k, k + 1)
                inside_chain &= (x_next ∈ E_by_k[k + 1])
            end

            if check_domain && hasproperty(run_result.concrete_system, :X)
                if !(x_next ∈ run_result.concrete_system.X)
                    left_domain = true
                    x = x_next
                    break
                end
            end

            x = x_next
        end

        final_in_target = (x ∈ E_target)
        success = inside_chain && final_in_target && !left_domain

        x_rollouts[s] = x_hist
        u_rollouts[s] = u_hist
        rollout_stats[s] = (;
            success = success,
            inside_chain = inside_chain,
            final_in_target = final_in_target,
            left_domain = left_domain,
        )
    end

    n_success = count(r -> r.success, rollout_stats)
    n_inside = count(r -> r.inside_chain, rollout_stats)
    n_target = count(r -> r.final_in_target, rollout_stats)
    n_left_domain = count(r -> r.left_domain, rollout_stats)

    summary = (;
        n_samples = n_samples,
        n_success = n_success,
        n_inside_chain = n_inside,
        n_final_in_target = n_target,
        n_left_domain = n_left_domain,
        success_rate = n_success / n_samples,
        inside_chain_rate = n_inside / n_samples,
        final_target_rate = n_target / n_samples,
        left_domain_rate = n_left_domain / n_samples,
    )

    if verbose
        println("\n=== Validation empirique des κ_k sur modele concret ===")
        println("samples=", summary.n_samples)
        println(
            "success=",
            summary.n_success,
            " (",
            round(100 * summary.success_rate; digits = 1),
            "%)",
        )
        println(
            "inside_chain=",
            summary.n_inside_chain,
            " (",
            round(100 * summary.inside_chain_rate; digits = 1),
            "%)",
        )
        println(
            "final_in_target=",
            summary.n_final_in_target,
            " (",
            round(100 * summary.final_target_rate; digits = 1),
            "%)",
        )
        println(
            "left_domain=",
            summary.n_left_domain,
            " (",
            round(100 * summary.left_domain_rate; digits = 1),
            "%)",
        )
    end

    return (;
        summary = summary,
        k_sequence = k_sorted,
        E_init = E_init,
        E_target = E_target,
        x0_samples = x0_samples,
        x_rollouts = x_rollouts,
        u_rollouts = u_rollouts,
        rollout_stats = rollout_stats,
    )
end
