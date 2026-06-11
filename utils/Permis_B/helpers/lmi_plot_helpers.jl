import LinearAlgebra as LA

"""
    print_state_list(state_list)

Affiche chaque etat de la trajectoire avec son indice.
"""
function print_state_list(state_list)
    println("State list (" * string(length(state_list)) * " states):")
    for (i, x) in enumerate(state_list)
        println("[" * string(i) * "] " * string(x))
    end
end

"""
    projeter_ellipsoide_2d(E; dims=(1, 2))

Projette une ellipsoide haute dimension sur un plan 2D via la matrice
de covariance `Q = inv(P)`, puis reconstruit la forme quadratique projetee.
"""
function projeter_ellipsoide_2d(E::UT.Ellipsoid; dims = (1, 2))
    i, j = dims
    P = Matrix(E.P)

    # E = {x | (x-c)'P(x-c) <= 1}, projection via la forme Q = inv(P).
    Q = try
        inv(P)
    catch
        LA.pinv(P)
    end
    Q2 = Q[[i, j], [i, j]]
    P2 = try
        inv(Q2)
    catch
        LA.pinv(Q2)
    end

    c2 = [E.c[i], E.c[j]]
    return UT.Ellipsoid(P2, c2)
end

"""
    plot_ellipsoide_terminale(state_list, E_terminal; output_dir=...)

Trace la trajectoire `(x1, x2)` et l'ellipsoide terminale projetee.
"""
function plot_ellipsoide_terminale(state_list, E_terminal; output_dir = dirname(@__DIR__))
    traj = ST.Trajectory(state_list)
    E_terminal_xy = projeter_ellipsoide_2d(E_terminal; dims = (1, 2))

    fig = plot(;
        aspect_ratio = :equal,
        legend = true,
        title = "Test LMI - Ellipsoide terminale",
    )
    plot!(
        fig,
        traj;
        dims = [1, 2],
        ms = 2.0,
        arrows = false,
        color = :blue,
        label = "Trajectoire",
    )
    plot!(
        fig,
        E_terminal_xy;
        color = :orange,
        opacity = 0.35,
        label = "E terminale proj. (x1,x2)",
    )
    scatter!(
        fig,
        [E_terminal.c[1]],
        [E_terminal.c[2]];
        color = :red,
        marker = :star5,
        ms = 7,
        label = "Centre E",
    )

    out_path = joinpath(output_dir, "lmi_ellipsoide_terminale_12.pdf")
    savefig(fig, out_path)
    display(fig)
    return out_path
end

"""
    plot_test_transition_backward(state_list, E_terminal, E_prev; output_dir=...)

Visualise une transition backward unique (avant-dernier vers terminal).
"""
function plot_test_transition_backward(
    state_list,
    E_terminal,
    E_prev;
    output_dir = dirname(@__DIR__),
)
    traj = ST.Trajectory(state_list)
    E_terminal_xy = projeter_ellipsoide_2d(E_terminal; dims = (1, 2))
    E_prev_xy = projeter_ellipsoide_2d(E_prev; dims = (1, 2))

    fig = plot(;
        aspect_ratio = :equal,
        legend = true,
        title = "Test LMI - Transition backward (avant-dernier -> terminal)",
    )
    plot!(
        fig,
        traj;
        dims = [1, 2],
        ms = 2.0,
        arrows = false,
        color = :blue,
        label = "Trajectoire",
    )
    plot!(fig, E_terminal_xy; color = :orange, opacity = 0.30, label = "E terminale")
    plot!(fig, E_prev_xy; color = :green, opacity = 0.30, label = "E avant-dernier")
    scatter!(
        fig,
        [E_terminal.c[1]],
        [E_terminal.c[2]];
        color = :red,
        marker = :star5,
        ms = 7,
        label = "Centre terminal",
    )
    scatter!(
        fig,
        [E_prev.c[1]],
        [E_prev.c[2]];
        color = :black,
        marker = :diamond,
        ms = 6,
        label = "Centre E_prev",
    )

    out_path = joinpath(output_dir, "lmi_transition_backward_12.pdf")
    savefig(fig, out_path)
    display(fig)
    return out_path
end

"""
    plot_transitions_backward_chaine(state_list, ellipsoides; output_dir=...)

Affiche toute la chaine backward d'ellipsoides dans un plan 2D choisi.
Par defaut `dims=(1,2)` (plan position). Pour les angles: `dims=(3,4)`.
"""
function plot_transitions_backward_chaine(
    state_list,
    ellipsoides;
    output_dir = dirname(@__DIR__),
    dims = (1, 2),
    title = nothing,
    filename = nothing,
)
    isempty(ellipsoides) && error("Aucune ellipsoide a afficher.")
    d1, d2 = dims

    traj = ST.Trajectory(state_list)
    nb = length(ellipsoides)
    cols = palette(:viridis, nb)
    title_txt =
        title === nothing ? "Test LMI - Chaine transitions backward ($(d1),$(d2))" : title

    fig = plot(; aspect_ratio = :equal, legend = true, title = title_txt)
    plot!(
        fig,
        traj;
        dims = [d1, d2],
        ms = 2.0,
        arrows = false,
        color = :blue,
        label = "Trajectoire",
    )

    for (i, E) in enumerate(ellipsoides)
        E_xy = projeter_ellipsoide_2d(E; dims = (d1, d2))
        lbl = i == 1 ? "E terminale" : (i == nb ? "E la plus reculee" : "")
        plot!(fig, E_xy; color = cols[i], opacity = 0.22, label = lbl)
    end

    xcent = [E.c[d1] for E in ellipsoides]
    ycent = [E.c[d2] for E in ellipsoides]
    scatter!(
        fig,
        xcent,
        ycent;
        color = :black,
        ms = 2.5,
        alpha = 0.8,
        label = "Centres ellipsoides",
    )

    filename_out =
        filename === nothing ? "lmi_transition_backward_$(d1)$(d2).pdf" : filename
    out_path = joinpath(output_dir, filename_out)
    savefig(fig, out_path)
    display(fig)
    return out_path
end

"""
    plot_kappa_rollouts_state_space(empirical_result; ...)

Trace les rollouts concrets echantillonnes avec les lois locales `κ_k`:
- trajectoires succes en vert,
- trajectoires echec en rouge,
- projection de l'ellipsoide initiale et finale.

Permet de verifier visuellement si les trajectoires echantillonnees
atteignent bien l'ensemble final.
"""
function plot_kappa_rollouts_state_space(
    empirical_result;
    output_dir = dirname(@__DIR__),
    dims = (1, 2),
    filename = nothing,
    title = nothing,
    domain = nothing,
    show_init::Bool = true,
    show_target::Bool = true,
)
    isempty(empirical_result.x_rollouts) && error("Aucune rollout a afficher.")
    d1, d2 = dims

    # Resume de succes (si disponible) pour enrichir le titre.
    succ_rate =
        hasproperty(empirical_result, :summary) &&
        hasproperty(empirical_result.summary, :success_rate) ?
        empirical_result.summary.success_rate : NaN
    title_txt =
        title !== nothing ? title :
        (
            isfinite(succ_rate) ?
            "Validation κ - rollouts concrets ($(d1),$(d2)) | succes=$(round(100*succ_rate; digits=1))%" :
            "Validation κ - rollouts concrets ($(d1),$(d2))"
        )

    fig = plot(; aspect_ratio = :equal, legend = true, title = title_txt)

    if domain !== nothing
        plot!(
            fig,
            domain;
            dims = [d1, d2],
            color = :grey,
            opacity = 0.12,
            label = "Domaine",
        )
    end

    if show_init && hasproperty(empirical_result, :E_init)
        E_init_xy = projeter_ellipsoide_2d(empirical_result.E_init; dims = (d1, d2))
        plot!(fig, E_init_xy; color = :deepskyblue3, opacity = 0.16, label = "E initiale")
    end

    if show_target && hasproperty(empirical_result, :E_target)
        E_target_xy = projeter_ellipsoide_2d(empirical_result.E_target; dims = (d1, d2))
        plot!(fig, E_target_xy; color = :orange, opacity = 0.30, label = "E finale")
        scatter!(
            fig,
            [empirical_result.E_target.c[d1]],
            [empirical_result.E_target.c[d2]];
            color = :orange,
            marker = :star5,
            ms = 6,
            label = "Centre E finale",
        )
    end

    # Points initiaux echantillonnes.
    if hasproperty(empirical_result, :x0_samples)
        xs0 = [x[d1] for x in empirical_result.x0_samples]
        ys0 = [x[d2] for x in empirical_result.x0_samples]
        scatter!(
            fig,
            xs0,
            ys0;
            color = :black,
            alpha = 0.35,
            marker = :xcross,
            ms = 2.0,
            label = "Echantillons initiaux",
        )
    end

    label_success_used = false
    label_fail_used = false
    for i in eachindex(empirical_result.x_rollouts)
        xseq = empirical_result.x_rollouts[i]
        traj = ST.Trajectory(xseq)
        is_success =
            hasproperty(empirical_result, :rollout_stats) &&
            i <= length(empirical_result.rollout_stats) ?
            empirical_result.rollout_stats[i].success : false

        col = is_success ? :seagreen3 : :crimson
        lbl =
            is_success ? (!label_success_used ? "Traj. succes" : "") :
            (!label_fail_used ? "Traj. echec" : "")
        plot!(
            fig,
            traj;
            dims = [d1, d2],
            ms = 1.0,
            arrows = false,
            color = col,
            alpha = 0.45,
            label = lbl,
        )
        if is_success
            label_success_used = true
        else
            label_fail_used = true
        end

        xend = xseq[end]
        scatter!(
            fig,
            [xend[d1]],
            [xend[d2]];
            color = col,
            alpha = 0.7,
            marker = :circle,
            ms = 2.0,
            label = "",
        )
    end

    filename_out =
        filename === nothing ? "kappa_rollouts_state_space_$(d1)$(d2).pdf" : filename
    out_path = joinpath(output_dir, filename_out)
    savefig(fig, out_path)
    display(fig)
    return out_path
end
