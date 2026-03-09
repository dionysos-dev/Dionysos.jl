using Plots
import StaticArrays: SVector
using Serialization
using SHA
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const OP = DI.Optim
const AB = OP.Abstraction
const SY = DI.Symbolic
using JuMP
import MathOptInterface as MOI
import IntervalArithmetic as IA
using MosekTools

#=
Permis_B.jl - Guide de lecture (vue d'ensemble)
================================================

Objectif
--------
Le script fait deux phases:
1) calcul d'une trajectoire nominale (abstraction + controleur),
2) diagnostic LMI backward sur cette trajectoire (mode symbolique uniquement).

Pipeline
--------
`main(...)` appelle:
- `run_nominal_simulation(...)`
- `run_lmi_diagnostic(...)`

Hypotheses
----------
- La simulation nominale utilise `AV.system(...)` (black-box).
- Les LMIs utilisent `AV.symbolic_system(...)` (discretisation configurable Euler/RK4).
- Les deux doivent rester coherents localement pour eviter les infeasibilites.

Reglage rapide
--------------
- Si echec LMI precoce: augmenter `rayon_terminal`, augmenter `λ`,
  et/ou reduire `ΔX`, `ΔU`.
=#

include(joinpath(@__DIR__, "helpers", "centered_system_helpers.jl"))
include(joinpath(@__DIR__, "helpers", "plot_helpers.jl"))
include(joinpath(@__DIR__, "helpers", "lmi_plot_helpers.jl"))
include(joinpath(@__DIR__, "helpers", "lmi_helpers.jl"))

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

Base.@kwdef struct SimulationConfig
    Δt::Float64 = 0.2
    hx::SVector{4, Float64} = SVector(0.4, 0.2, 5 * (pi / 180), 3 * (pi / 180))
    periodic_dims::SVector{2, Int} = SVector(3, 4)
    periodic_periods::SVector{2, Float64} = SVector(2pi, 2pi)
    periodic_start::SVector{2, Float64} = SVector(-pi, -pi)
    nstep::Int = 700
end

Base.@kwdef struct LMIConfig
    rayon_terminal::Float64 = 0.45
    λ::Float64 = 0.001
    maxδx::Float64 = 100.0
    maxδu::Float64 = 200.0
    # Optional external SDP optimizer override (e.g. MosekTools.Optimizer).
    # If `nothing`, Clarabel + clarabel_* settings are used.
    sdp_opt::Any = build_mosek_optimizer(verbose=false)
    # Nombre de sous-pas RK4 pour le modele symbolique.
    symbolic_rk4_substeps::Int = 1
    # LMI model switches:
    # - use_noise = false ignores additive disturbance term ν[j]
    # - use_logdet = false removes LogDetCone (no exponential cones)
    # - trace_reward enables objective λJ - α tr(L) when use_logdet=false
    ΔX = IA.IntervalBox(
        IA.interval(-0.1, 0.1),  # x1
        IA.interval(-0.1, 0.1),  # x2
        IA.interval(-0.1, 0.1),   # θ
        IA.interval(-0.1, 0.1)   # ϕ
    )
    ΔU = IA.IntervalBox(
        IA.interval(-0.2, 0.2),# v
        IA.interval(-0.2, 0.2) # δ (steering) mainteannt c'est tan(δ) (ça devrait etre largement)
    )
    ΔW::Any = IA.IntervalBox(IA.interval(0.0, 0.0), 1)
    verbose::Bool = true
end

function build_mosek_optimizer(; verbose::Bool = true)
    return optimizer_with_attributes(
        MosekTools.Optimizer,
        MOI.Silent() => !verbose,
        # 10 = detailed log, 0 = silent.
        MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => (verbose ? 10 : 0),
    )
end

function build_concrete_system()
    x_domain = UT.HyperRectangle( # backward + forward
        SVector(-1.0, -1.0, -pi, -pi),
        SVector(10.0, 9.0, pi, pi),
    )
    
    x_domain = AV.with_phi_limit(x_domain; phi_max = deg2rad(50.0))
    obstacles_xy = [ # c'est bout le benchmark classique
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 6.0), SVector(10.0, 9.0)),
    ]
    # obstacles_xy = [
    #     UT.HyperRectangle(SVector(-1.0, -1.0), SVector(10.0, 0.8)),  # trottoir/bordure
    #     UT.HyperRectangle(SVector(5.2, 0.8), SVector(10.0, 2.6)),     # voiture derrière
    #     UT.HyperRectangle(SVector(-2.0, 0.8), SVector(0.5, 2.6)),    # voiture devant
    # ]

    x_domain = AV.with_xy_obstacles(x_domain; obstacles2d = obstacles_xy)

    δ_max = pi/4
    σ_max = tan(δ_max)
    u_domain = UT.HyperRectangle(SVector(-5.0, -σ_max), SVector(5.0, σ_max)) # mixte

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(x_domain; _U_ = u_domain, params = params)
    return (; x_domain, u_domain, params, concrete_system)
end

function build_input_domain()
    inputs_delta = [
        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
    ]
    inputs = [[u[1], tan(u[2])] for u in inputs_delta]
    return DO.CustomList(inputs)
end

function build_control_problem()
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    initial_set = UT.HyperRectangle(
        SVector(-1.0, -1.0, -0.4, -0.4),
        SVector(1.0, 1.0, 0.4, 0.4),
    )
    target_set = UT.HyperRectangle( # il va faire de la marche avant
        SVector(9.0, 5.0, -5 * (pi / 180), -5 * (pi / 180)),
        SVector(10.0, 6.0, 5 * (pi / 180), 5 * (pi / 180)),
    )
    target_set = UT.HyperRectangle( # il va faire de la marche arrière
        SVector(9.0, 5.0, pi - 5 * (pi / 180), -5 * (pi / 180)),
        SVector(10.0, 6.0, pi + 5 * (pi / 180),  5 * (pi / 180)),
    )

    # pour le créneau arrière
    # x0 = SVector(9.2, 4.2, 0.0, 0.0)
    # initial_set = UT.HyperRectangle(
    #     SVector(8.8, 3.8, pi - 8*(pi/180), - 8*(pi/180)),
    #     SVector(9.6, 4.6, pi + 8*(pi/180), pi + 8*(pi/180)),
    # )
    # target_set = UT.HyperRectangle(
    #     SVector(2.2, 1.2, -5 * (pi / 180), -5 * (pi / 180)),
    #     SVector(4.0, 2.2,  5 * (pi / 180), 5 * (pi / 180)),
    # )


    return (; x0, initial_set, target_set)
end

periodicity_kwargs(cfg::SimulationConfig) = (;
    with_period = true,
    periodic_dims = cfg.periodic_dims,
    periodic_periods = cfg.periodic_periods,
    periodic_start = cfg.periodic_start,
)

"""
    unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)

Construit une copie "depliee" (unwrap) de `state_list` sur les dimensions
periodiques. Pour chaque angle, on remplace l'increment brut Δ par son
representant de plus petite norme:
`Δ_min = Δ - round(Δ / p) * p`, avec `p` la periode.

Mathematiquement, cela impose une trajectoire continue sur le revetement
universel de S¹, ce qui evite les faux sauts de ±2π dans les LMIs locales.
"""
function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims et periodic_periods doivent avoir la meme longueur.")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    println(periodic_dims)

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        1 <= d <= nx || error("Dimension periodique invalide: $d.")
        p > 0 || error("Periode invalide (<= 0): $p.")

        for k in 2:length(xs)
            Δ_raw = xs[k][d] - xs[k - 1][d]
            Δ_min = Δ_raw - round(Δ_raw / p) * p
            xs[k][d] = xs[k - 1][d] + Δ_min
        end
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

"""
    prepare_run_result_for_lmi(run_result; ...)

Retourne une vue `run_result` adaptee a la phase LMI:
- soit inchangée,
- soit avec `state_list` unwrap sur les dimensions periodiques.
"""
function prepare_run_result_for_lmi(
    run_result;
    enable_unwrap::Bool = true,
    periodic_dims = nothing,
    periodic_periods = nothing,
)
    enable_unwrap || return run_result
    isempty(run_result.state_list) && return run_result

    dims =
        periodic_dims === nothing ?
        (hasproperty(run_result, :periodic_dims) ? getproperty(run_result, :periodic_dims) : nothing) :
        periodic_dims
    periods =
        periodic_periods === nothing ?
        (hasproperty(run_result, :periodic_periods) ? getproperty(run_result, :periodic_periods) :
         nothing) : periodic_periods

    (dims === nothing || periods === nothing || isempty(dims)) && return run_result

    unwrapped_states = unwrap_periodic_state_list(run_result.state_list, dims, periods)
    return (; run_result..., state_list = unwrapped_states)
end

function nominal_cache_key(cfg::SimulationConfig)
    payload = (
        Δt = cfg.Δt,
        hx = Tuple(cfg.hx),
        periodic_dims = Tuple(cfg.periodic_dims),
        periodic_periods = Tuple(cfg.periodic_periods),
        periodic_start = Tuple(cfg.periodic_start),
        nstep = cfg.nstep,
    )
    return bytes2hex(sha1(repr(payload)))
end

function nominal_cache_path(
    cfg::SimulationConfig;
    cache_dir::AbstractString = joinpath(@__DIR__, "cache_nominal"),
)
    mkpath(cache_dir)
    return joinpath(cache_dir, "nominal_" * nominal_cache_key(cfg) * ".jls")
end

function run_nominal_simulation(
    cfg::SimulationConfig = SimulationConfig();
    save_plots::Bool = false,
    show_animation::Bool = false,
    centered_animation_gif::Union{Nothing, AbstractString} = nothing,
    centered_animation_title::AbstractString = "centered_simu",
    use_cache::Bool = true,
    force_recompute::Bool = false,
    cache_dir::AbstractString = joinpath(@__DIR__, "cache_nominal"),
)
    system_cfg = build_concrete_system()
    centered_gif_path =
        centered_animation_gif === nothing ?
        joinpath(@__DIR__, "centered_simu.gif") : String(centered_animation_gif)
    cache_file = nominal_cache_path(cfg; cache_dir = cache_dir)
    if use_cache && !force_recompute && isfile(cache_file)
        try
            payload = deserialize(cache_file)
            state_list = payload.state_list
            input_list = payload.input_list
            println("Trajectoire nominale chargee depuis le cache: ", cache_file)
            animation = nothing
            if show_animation
                x_traj = ST.Trajectory(state_list)
                u_traj = ST.Trajectory(input_list)
                animation = plot_articulated_vehicle!(
                    AV,
                    system_cfg.concrete_system,
                    system_cfg.params,
                    x_traj,
                    u_traj;
                    every = 1,
                    dt = 0.09,
                    giffile = centered_gif_path,
                    title = centered_animation_title,
                )
                println("Animation centered sauvegardee: ", centered_gif_path)
            end
            return (;
                state_list,
                input_list,
                animation,
                concrete_system = system_cfg.concrete_system,
                Δt = cfg.Δt,
                params = system_cfg.params,
                periodic_dims = cfg.periodic_dims,
                periodic_periods = cfg.periodic_periods,
                periodic_start = cfg.periodic_start,
                cache_file,
                centered_animation_gif = show_animation ? centered_gif_path : nothing,
                from_cache = true,
            )
        catch err
            println("Cache nominal illisible, recalcul force (", err, ").")
        end
    end

    control_cfg = build_control_problem()
    Udom = build_input_domain()

    opt = build_uniform_grid_abstraction(
        system_cfg.concrete_system,
        cfg.Δt,
        cfg.hx,
        Udom,
        AV.jacobian_bound(system_cfg.params);
        periodicity_kwargs(cfg)...,
    )

    controller, target_set = build_uniform_grid_controller!(
        opt,
        system_cfg.concrete_system,
        control_cfg.initial_set,
        control_cfg.target_set,
    )

    x_traj, u_traj = simulate_closed_loop(
        system_cfg.concrete_system,
        controller,
        cfg.Δt,
        control_cfg.x0,
        target_set;
        periodicity_kwargs(cfg)...,
        nstep = cfg.nstep,
    )
    state_list, input_list = collect_trajectory_lists(x_traj, u_traj)

    if save_plots
        save_state_space_plots!(
            @__DIR__,
            opt,
            system_cfg.concrete_system,
            control_cfg.initial_set,
            control_cfg.target_set,
            x_traj;
            periodicity_kwargs(cfg)...,
        )
    end

    animation =
        show_animation ? plot_articulated_vehicle!(
            AV,
            system_cfg.concrete_system,
            system_cfg.params,
            x_traj,
            u_traj;
            every = 1,
            dt = 0.09,
            giffile = centered_gif_path,
            title = centered_animation_title,
        ) : nothing
    if show_animation
        println("Animation centered sauvegardee: ", centered_gif_path)
    end

    if use_cache
        serialize(cache_file, (; state_list, input_list))
        println("Trajectoire nominale sauvegardee dans le cache: ", cache_file)
    end

    return (;
        state_list,
        input_list,
        animation,
        concrete_system = system_cfg.concrete_system,
        Δt = cfg.Δt,
        params = system_cfg.params,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        cache_file,
        centered_animation_gif = show_animation ? centered_gif_path : nothing,
        from_cache = false,
    )
end

function build_symbolic_lmi_context(run_result, cfg::LMIConfig)
    sym_sys = AV.symbolic_system(
        run_result.concrete_system.X;
        _U_ = run_result.concrete_system.U,
        params = run_result.params,
        Ts = run_result.Δt,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rk4_num_substeps = cfg.symbolic_rk4_substeps,
    )

    nx = length(run_result.state_list[1])
    nu = isempty(run_result.input_list) ? UT.get_dims(run_result.concrete_system.U) :
         length(run_result.input_list[1])
    S = matrice_cout_identite(nx, nu)

    return construire_contexte_lmi(;
        fsymbolic = sym_sys.fsymbolic,
        x = sym_sys.x,
        u = sym_sys.u,
        w = sym_sys.w,
        ΔX = sym_sys.ΔX,
        ΔU = sym_sys.ΔU,
        ΔW = sym_sys.ΔW,
        Uformat = sym_sys.Uformat,
        Wformat = sym_sys.Wformat,
        transition_cost = S,
        sdp_opt = cfg.sdp_opt,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        λ = cfg.λ
    )
end

function run_lmi_diagnostic(
    run_result;
    cfg::LMIConfig = LMIConfig(),
    output_dir = @__DIR__,
    save_plots::Bool = false,
    unwrap_periodic_angles::Bool = true,
    periodic_dims = nothing,
    periodic_periods = nothing,
)
    run_result_lmi = prepare_run_result_for_lmi(
        run_result;
        enable_unwrap = unwrap_periodic_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
    )

    println("\n=== Test LMI (mode symbolique) ===")
    E_terminal = creer_ellipsoide_terminal(run_result_lmi.state_list; rayon = cfg.rayon_terminal)
    println("Ellipsoide terminale creee.")
    println("Centre: ", E_terminal.c)
    println("P[1,1]: ", E_terminal.P[1, 1])

    pdf_terminal = nothing
    if save_plots
        pdf_terminal = plot_ellipsoide_terminale(
            run_result_lmi.state_list,
            E_terminal;
            output_dir = output_dir,
        )
        println("Plot enregistre: ", pdf_terminal)
    end

    ctx = build_symbolic_lmi_context(run_result_lmi, cfg)
    transitions = synthetiser_transitions_backward(
        run_result_lmi.state_list,
        run_result_lmi.input_list,
        ctx;
        E_terminal = E_terminal,
        verbose = cfg.verbose,
    )

    println("Transitions tentees: ", length(run_result_lmi.state_list) - 1)
    println("Transitions reussies: ", length(transitions.indices))
    println("Ellipsoides calculees (avec terminale): ", length(transitions.ellipsoides))
    if transitions.success
        println("Chaine backward complete reussie.")
    else
        println("Chaine backward incomplete. premier echec k = ", transitions.failed_k)
    end
    if !isempty(transitions.couts)
        println("Cout min/max = ", minimum(transitions.couts), " / ", maximum(transitions.couts))
    end

    pdf_transitions = nothing
    pdf_transitions_angles = nothing
    if save_plots
        pdf_transitions = plot_transitions_backward_chaine(
            run_result_lmi.state_list,
            transitions.ellipsoides;
            output_dir = output_dir,
            dims = (1, 2),
        )
        println("Plot transition enregistre: ", pdf_transitions)

        pdf_transitions_angles = plot_transitions_backward_chaine(
            run_result_lmi.state_list,
            transitions.ellipsoides;
            output_dir = output_dir,
            dims = (3, 4),
            title = "Test LMI - Chaine transitions backward (angles theta, phi)",
            filename = "lmi_transition_backward_34_angles.pdf",
        )
        println("Plot transition angles enregistre: ", pdf_transitions_angles)
    end

    return (; E_terminal, transitions, pdf_terminal, pdf_transitions, pdf_transitions_angles)
end

"""
    run_kappa_empirical_check(run_result, diag_result; cfg, n_samples, ...)

Wrapper pratique pour valider empiriquement les lois `κ_k` synthetisees
par `run_lmi_diagnostic` sur le modele concret.
"""
function run_kappa_empirical_check(
    run_result,
    diag_result;
    cfg::LMIConfig = LMIConfig(),
    n_samples::Int = 100,
    num_substeps::Int = 1,
    project_inputs::Bool = true,
    check_domain::Bool = true,
    verbose::Bool = true,
)
    ctx = build_symbolic_lmi_context(run_result, cfg)
    return simuler_kappa_sur_modele_concret(
        run_result,
        diag_result.transitions,
        ctx;
        n_samples = n_samples,
        num_substeps = num_substeps,
        project_inputs = project_inputs,
        check_domain = check_domain,
        verbose = verbose,
    )
end

"""
    animate_kappa_rollout(run_result, empirical_result; sample_idx, ...)

Anime une rollout issue de `run_kappa_empirical_check` en reutilisant
le meme rendu que la centered simulation (`plot_articulated_vehicle!`).
"""
function animate_kappa_rollout(
    run_result,
    empirical_result;
    sample_idx::Int = 1,
    giffile::Union{Nothing, String} = nothing,
    fps::Int = 20,
    every::Int = 1,
    dt::Float64 = 0.09,
    title::String = "ellipsoidal controller",
)
    n = length(empirical_result.x_rollouts)
    (1 <= sample_idx <= n) ||
        error("sample_idx doit etre dans 1:$n.")

    x_seq = empirical_result.x_rollouts[sample_idx]
    u_seq = empirical_result.u_rollouts[sample_idx]

    x_traj = ST.Trajectory(x_seq)
    u_traj = ST.Trajectory(u_seq)

    return plot_articulated_vehicle!(
        AV,
        run_result.concrete_system,
        run_result.params,
        x_traj,
        u_traj;
        giffile = giffile,
        fps = fps,
        every = every,
        dt = dt,
        title = title,
    )
end

"""
    plot_kappa_empirical_state_space(run_result, empirical_result; ...)

Wrapper pratique pour tracer les trajectoires concretes issues de
`run_kappa_empirical_check` dans l'espace d'etat.
"""
function plot_kappa_empirical_state_space(
    run_result,
    empirical_result;
    output_dir = @__DIR__,
    dims = (1, 2),
    filename = nothing,
    title = nothing,
)
    return plot_kappa_rollouts_state_space(
        empirical_result;
        output_dir = output_dir,
        dims = dims,
        filename = filename,
        title = title,
        domain = run_result.concrete_system.X,
    )
end

function script(;
    save_plots = true,
    show_animation = false,
    use_cache = true,
    force_recompute = false,
)
    return run_nominal_simulation(;
        save_plots = save_plots,
        show_animation = show_animation,
        use_cache = use_cache,
        force_recompute = force_recompute,
    )
end

function main(;
    save_plots = true,
    show_animation = true,
    run_lmi = true,
    run_lmi_empirical_check::Bool = true,
    lmi_check_n_samples::Int = 100,
    lmi_check_num_substeps::Int = 5,
    lmi_check_project_inputs::Bool = true,
    lmi_check_check_domain::Bool = true,
    lmi_check_verbose::Bool = true,
    show_lmi_animation::Bool = true,
    lmi_animation_sample_idx::Int = 1,
    lmi_animation_n_samples::Int = 20,
    lmi_animation_gif::Union{Nothing, String} = nothing,
    lmi_animation_title::String = "ellipsoidal controller",
    use_cached_nominal = false,
    force_recompute_nominal = false,
    nominal_cache_dir = joinpath(@__DIR__, "cache_nominal"),
    sim_cfg::SimulationConfig = SimulationConfig(),
    lmi_cfg::LMIConfig = LMIConfig(),
)
    run_result = run_nominal_simulation(
        sim_cfg;
        save_plots = save_plots,
        show_animation = show_animation,
        use_cache = use_cached_nominal,
        force_recompute = force_recompute_nominal,
        cache_dir = nominal_cache_dir,
    )
    print_state_list(run_result.state_list)

    lmi_result = nothing
    lmi_empirical_result = nothing
    lmi_animation = nothing
    if run_lmi
        try
            lmi_result = run_lmi_diagnostic(
                run_result;
                cfg = lmi_cfg,
                output_dir = @__DIR__,
                save_plots = save_plots,
            )
        catch err
            println("Echec phase LMI: ", err)
        end
    end

    # Test concret des lois κ_k (closed-loop sur modele discretise du systeme concret).
    # Ce test peut etre active seul (sans animation), ou implicitement via
    # show_lmi_animation=true.
    do_empirical_check = run_lmi_empirical_check || show_lmi_animation
    if do_empirical_check
        if lmi_result === nothing
            println("Test concret κ ignore: aucun resultat LMI disponible.")
        elseif isempty(lmi_result.transitions.indices)
            println("Test concret κ ignore: aucune transition backward valide.")
        else
            n_samples_for_check =
                show_lmi_animation ?
                max(max(lmi_animation_n_samples, lmi_animation_sample_idx), lmi_check_n_samples) :
                lmi_check_n_samples
            lmi_empirical_result = run_kappa_empirical_check(
                run_result,
                lmi_result;
                cfg = lmi_cfg,
                n_samples = n_samples_for_check,
                num_substeps = lmi_check_num_substeps,
                project_inputs = lmi_check_project_inputs,
                check_domain = lmi_check_check_domain,
                verbose = lmi_check_verbose,
            )

            if show_lmi_animation
                lmi_gif_path =
                    lmi_animation_gif === nothing ?
                    joinpath(@__DIR__, "ellipsoidal_controller.gif") : String(lmi_animation_gif)
                lmi_animation = animate_kappa_rollout(
                    run_result,
                    lmi_empirical_result;
                    sample_idx = lmi_animation_sample_idx,
                    giffile = lmi_gif_path,
                    title = lmi_animation_title,
                )
                println("Animation LMI sauvegardee: ", lmi_gif_path)
            end
        end
    end

    return (; run_result, lmi_result, lmi_empirical_result, lmi_animation)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
