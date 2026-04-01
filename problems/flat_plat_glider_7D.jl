module FlatPlateGlider7D

# Benchmark flat-plate glider 7D.
# Source principale : Moore12.pdf, section IV "AIRCRAFT MODEL".

using StaticArrays
using MathematicalSystems
using Dionysos
using Plots
import Symbolics
import IntervalArithmetic as IA

const UT = Dionysos.Utils
const ST = Dionysos.System
const PB = Dionysos.Problem
const SY = Dionysos.Symbolic

# ----------------------------
# Parametres
# ----------------------------
Base.@kwdef struct Params{T}
    m::T = 0.085                  # masse (kg). Valeur explicite : 85 g dans la description experimentale.
    I::T = 5.0e-4                # inertie de tangage (kg.m^2). Non donnee explicitement dans la section IV accessible : valeur benchmark a verifier.
    g::T = 9.81                  # gravite (m/s^2).
    rho::T = 1.225               # densite de l'air (kg/m^3) au niveau de la mer.
    Sw::T = 0.085                # surface aerodynamique de l'aile (m^2). La repartition exacte aile/elevator n'est pas explicite dans l'article : a ajuster si necessaire.
    Se::T = 0.012                # surface aerodynamique de l'elevator (m^2). Choix benchmark compatible avec une surface totale proche de 0.102 m^2.
    lw::T = 0.03                 # distance COM -> centre aero aile (m). Non explicite dans l'article : valeur benchmark a verifier.
    l::T = 0.23                  # distance COM -> charniere elevator (m). Non explicite dans l'article : valeur benchmark a verifier.
    le::T = 0.05                 # distance charniere -> centre aero elevator (m). Non explicite dans l'article : valeur benchmark a verifier.

    phi_min::T = -pi / 3         # borne article / benchmark sur l'angle elevator.
    phi_max::T = pi / 8          # borne article / benchmark sur l'angle elevator.
    phi_dot_max::T = 13.0        # borne article / benchmark sur la vitesse de l'elevator.
    vx_abs_bound::T = 10.0       # borne benchmark pour |vx|, utilisee par les boites finies et les majorations conservatives.
    vz_abs_bound::T = 8.0        # borne benchmark pour |vz|, utilisee par les boites finies et les majorations conservatives.
    q_abs_bound::T = 20.0        # borne benchmark pour |q|, utilisee par les boites finies et les majorations conservatives.

    jacobian_fd_step::T = 1.0e-6      # pas de differences finies pour la jacobienne locale. Parametre purement numerique, expose pour etre ajustable.
    bound_velocity_probe::T = 0.25    # petite vitesse de sondage pour les bornes globales, afin d'eviter d'echantillonner exactement la singularite de atan(0,0).
    jacobian_bound_inflation::T = 1.25 # facteur de securite sur la borne jacobienne issue de l'echantillonnage.
end

# ----------------------------
# Dynamique continue : xdot = f(x,u)
# x = [px, pz, θ, φ, vx, vz, q]
# u = [φdot]
# ----------------------------

# Petit helper local :
# on factorise ici les quantites aerodynamiques communes a la dynamique continue
# et au systeme symbolique. Cela garde le code lisible sans construire un framework.
function _aero_terms(θ, φ, vx, vz, q, φdot, p::Params)
    θe = θ + φ

    # Normales aux surfaces, exactement dans l'esprit de la section IV.
    nwx = -sin(θ)
    nwz = cos(θ)
    nex = -sin(θe)
    nez = cos(θe)

    # Vitesses des centres aerodynamiques.
    vxw = vx + p.lw * q * sin(θ)
    vzw = vz - p.lw * q * cos(θ)

    vxe = vx + p.l * q * sin(θ) + p.le * (q + φdot) * sin(θe)
    vze = vz - p.l * q * cos(θ) - p.le * (q + φdot) * cos(θe)

    # atan(y, x) conserve correctement le quadrant, ce qui est important ici.
    αw = θ - atan(vzw, vxw)
    αe = θe - atan(vze, vxe)

    # Le cahier des charges demande ici la forme scalaire simplifiee suivante.
    fw = p.rho * p.Sw * (vxw^2 + vzw^2) * sin(αw)
    fe = p.rho * p.Se * (vxe^2 + vze^2) * sin(αe)

    return (; θe, nwx, nwz, nex, nez, vxw, vzw, vxe, vze, αw, αe, fw, fe)
end

function dynamic(p::Params = Params())
    return (x, u) -> begin
        px, pz, θ, φ, vx, vz, q = x
        φdot = u[1]

        aero = _aero_terms(θ, φ, vx, vz, q, φdot, p)

        # Remarque importante :
        # le PDF OCR est ambigu sur le terme "lcφ".
        # Ici on suit l'interpretation demandee : l*cos(φ) + le.
        qdot = (-aero.fw * p.lw - aero.fe * (p.l * cos(φ) + p.le)) / p.I

        return SVector{7}(
            vx,
            vz,
            q,
            φdot,
            (aero.fw * aero.nwx + aero.fe * aero.nex) / p.m,
            (aero.fw * aero.nwz + aero.fe * aero.nez - p.m * p.g) / p.m,
            qdot,
        )
    end
end

# ----------------------------
# Jacobienne A(x,u) = ∂f/∂x
# ----------------------------

# On privilegie ici une routine locale de differences centrees :
# l'expression analytique exacte est longue, fragile et peu relisible pour ce benchmark.
function _numeric_jacobian_x(f, x, u, p::Params)
    xv = SVector{7, Float64}(x)
    uv = SVector{1, Float64}(u)
    h0 = Float64(p.jacobian_fd_step)

    A = Matrix{Float64}(undef, 7, 7)
    for i in 1:7
        hi = max(h0, h0 * abs(xv[i]))
        dx = SVector{7, Float64}(ntuple(k -> k == i ? hi : 0.0, 7))
        fp = f(xv + dx, uv)
        fm = f(xv - dx, uv)
        A[:, i] .= (fp - fm) / (2 * hi)
    end

    return SMatrix{7, 7, Float64, 49}(A)
end

function jacobian(p::Params = Params())
    f = dynamic(p)
    return (x, u) -> _numeric_jacobian_x(f, x, u, p)
end

# ----------------------------
# Bornes jacobiennes / Lipschitz
# ----------------------------

function _phi_samples(p::Params)
    ϕmid = 0.5 * (p.phi_min + p.phi_max)
    return (Float64(p.phi_min), Float64(ϕmid), Float64(p.phi_max))
end

function _signed_probe_samples(maxabs, probe)
    a = Float64(abs(maxabs))
    if iszero(a)
        return (0.0,)
    end

    b = min(Float64(abs(probe)), a)
    vals = sort(unique([-a, -b, b, a]))
    return Tuple(vals)
end

function _centered_samples(maxabs)
    a = Float64(abs(maxabs))
    if iszero(a)
        return (0.0,)
    end
    return (-a, 0.0, a)
end

function _global_jacobian_bound_matrix(p::Params)
    jac = jacobian(p)

    # Borne conservative :
    # on balaie une petite grille de benchmark sur les variables qui influencent
    # vraiment la dynamique. px et pz n'apparaissent pas explicitement dans f.
    θ_samples = (-pi, -pi / 2, 0.0, pi / 2, pi)
    φ_samples = _phi_samples(p)
    vx_samples = _signed_probe_samples(p.vx_abs_bound, p.bound_velocity_probe)
    vz_samples = _signed_probe_samples(p.vz_abs_bound, p.bound_velocity_probe)
    q_samples = _centered_samples(p.q_abs_bound)
    φdot_samples = (-Float64(p.phi_dot_max), 0.0, Float64(p.phi_dot_max))

    Amax = zeros(7, 7)
    for θ in θ_samples
        for φ in φ_samples
            for vx in vx_samples
                for vz in vz_samples
                    for q in q_samples
                        for φdot in φdot_samples
                            x = SVector(0.0, 0.0, θ, φ, vx, vz, q)
                            u = SVector(φdot)
                            A = jac(x, u)
                            Amax .= max.(Amax, abs.(Matrix(A)))
                        end
                    end
                end
            end
        end
    end

    Amax .*= Float64(p.jacobian_bound_inflation)
    return SMatrix{7, 7, Float64, 49}(Amax)
end

function jacobian_bound(p::Params = Params())
    Jb = _global_jacobian_bound_matrix(p)
    return _u -> Jb
end

function bound_norm_jacobian(p::Params = Params())
    Jb = _global_jacobian_bound_matrix(p)
    # Norme infinie matricielle : monotone et compatible avec une lecture "croissance coordonnee par coordonnee".
    Lj = maximum(vec(sum(Matrix(Jb), dims = 2)))
    return _u -> Lj
end

function bound_norm_hessian_tensor(p::Params = Params())
    Jb = _global_jacobian_bound_matrix(p)
    Lj = maximum(vec(sum(Matrix(Jb), dims = 2)))

    # Borne tres conservative :
    # on ne cherche pas une borne fine du tenseur hessien, seulement un majorant
    # monotone et exploitable pour des heuristiques locales d'abstraction.
    Vbody = hypot(p.vx_abs_bound, p.vz_abs_bound)
    Vwing = Vbody + abs(p.lw) * p.q_abs_bound
    Vtail = Vbody + abs(p.l) * p.q_abs_bound + abs(p.le) * (p.q_abs_bound + p.phi_dot_max)
    Vscale = max(Vwing, Vtail, p.bound_velocity_probe)

    # Echelle grossiere du gain aerodynamique sur les lignes acceleration / moment.
    accel_gain = p.rho * (p.Sw + p.Se) * Vscale / p.m
    lever_max = max(
        abs(p.lw),
        abs(p.l * cos(p.phi_min) + p.le),
        abs(p.l * cos(p.phi_max) + p.le),
    )
    moment_gain = p.rho * (p.Sw + p.Se) * Vscale * max(lever_max, abs(p.le)) / p.I

    Hbound = Lj * (accel_gain + moment_gain)
    return _u -> Hbound
end

# ----------------------------
# Wrapper systeme Dionysos
# ----------------------------
function system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(-13.0), SVector(13.0)),
    params::Params = Params(),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

"""
    symbolic_system(_X_; _U_, params, Ts, ΔX, ΔU, ΔW)

Version symbolique discrete pour la synthese LMI.
Le systeme black-box reste utilise pour la simulation nominale.
Le systeme symbolique sert ici a la linearisation et aux outils LMI.

"""
function symbolic_system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(-13.0), SVector(13.0)),
    params::Params = Params(),
    Ts::Float64 = 0.01,
    ΔX = IA.IntervalBox(IA.interval(-0.2, 0.2), 7),
    ΔU = IA.IntervalBox(IA.interval(-0.15, 0.15), 1),
    ΔW = IA.IntervalBox(IA.interval(0.0, 0.0), 1),
    rk4_num_substeps::Int = 1,
    obstacles = Any[],
)
    Symbolics.@variables px pz θ φ vx vz q φdot w1 T
    x = [px; pz; θ; φ; vx; vz; q]
    u = [φdot]
    w = [w1]

    function f_cont_expr(xloc, uloc)
        θloc = xloc[3]
        φloc = xloc[4]
        vxloc = xloc[5]
        vzloc = xloc[6]
        qloc = xloc[7]
        φdotloc = uloc[1]

        aero = _aero_terms(θloc, φloc, vxloc, vzloc, qloc, φdotloc, params)

        # Meme remarque qu'en dynamique numerique :
        # le terme OCR "lcφ" est interprete comme l*cos(φ) + le.
        qdot = (-aero.fw * params.lw - aero.fe * (params.l * cos(φloc) + params.le)) / params.I

        return [
            vxloc
            vzloc
            qloc
            φdotloc
            (aero.fw * aero.nwx + aero.fe * aero.nex) / params.m
            (aero.fw * aero.nwz + aero.fe * aero.nez - params.m * params.g) / params.m
            qdot
        ]
    end

    # Discretisation symbolique choisie pour la linearisation LMI.
    f_disc = ST.runge_kutta4(f_cont_expr, x, u, T, rk4_num_substeps)

    fsymbolicT = eval(ST.build_function(f_disc, x, u, w, T)[1])
    fsymbolic = Symbolics.substitute(f_disc, Dict(T => Ts))

    # Pas de bruit additif dans ce benchmark, mais l'API LMI de Dionysos
    # attend tout de meme un format de bruit.
    Wset = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    Uformat = UT.format_input_set(_U_)
    Wformat = UT.format_noise_set(Wset)

    f_cont_fun = dynamic(params)
    function f_eval(xv, uv, _wv)
        xsv = SVector{7, Float64}(xv)
        usv = SVector{1, Float64}(uv)
        xnext = ST.runge_kutta4(f_cont_fun, xsv, usv, Ts, rk4_num_substeps)
        return collect(xnext)
    end
    function f_backward_eval(xv, uv, _wv)
        xsv = SVector{7, Float64}(xv)
        usv = SVector{1, Float64}(uv)
        xprev = ST.runge_kutta4(f_cont_fun, xsv, usv, -Ts, rk4_num_substeps)
        return collect(xprev)
    end

    return ST.SymbolicSystem(
        fsymbolicT,
        fsymbolic,
        Ts,
        length(x),
        length(u),
        length(w),
        x,
        u,
        w,
        ΔX,
        ΔU,
        ΔW,
        _X_,
        _U_,
        Wset,
        obstacles,
        f_eval,
        f_backward_eval,
        Uformat,
        Wformat,
    )
end

# ----------------------------
# Helpers de domaine
# ----------------------------
function with_phi_limit(_X_::UT.HyperRectangle; phi_min = -pi / 3, phi_max = pi / 8)
    lb = SVector(_X_.lb[1], _X_.lb[2], _X_.lb[3], phi_min, _X_.lb[5], _X_.lb[6], _X_.lb[7])
    ub = SVector(_X_.ub[1], _X_.ub[2], _X_.ub[3], phi_max, _X_.ub[5], _X_.ub[6], _X_.ub[7])
    return UT.HyperRectangle(lb, ub)
end

# ----------------------------
# Factory de probleme benchmark
# ----------------------------
function problem(;
    params::Params = Params(),
    transition_cost = nothing,
    state_cost = nothing,
    terminal_cost = PB.Infinity(),
)
    # Domaine benchmark fini :
    # ces bornes ne pretendent pas etre "la" boite exacte de l'article.
    # Elles sont choisies pour contenir la manoeuvre de perch avec une marge
    # tout en restant compatibles avec les HyperRectangles de Dionysos.
    _X_ = UT.HyperRectangle(
        SVector(
            -4.5,
            -1.5,
            -pi / 2,
            params.phi_min,
            -params.vx_abs_bound,
            -params.vz_abs_bound,
            -params.q_abs_bound,
        ),
        SVector(
            0.5,
            1.5,
            pi,
            params.phi_max,
            params.vx_abs_bound,
            params.vz_abs_bound,
            params.q_abs_bound,
        ),
    )

    # Petite boite autour de l'etat initial nominal de l'article.
    x0 = SVector(3.5, 0.1, 0.0, 0.0, 7.0, 0.0, 0.0)
    _I_ = UT.HyperRectangle(
        x0 - SVector(0.05, 0.05, 0.05, 0.02, 0.2, 0.2, 0.5),
        x0 + SVector(0.05, 0.05, 0.05, 0.02, 0.2, 0.2, 0.5),
    )

    # Region terminale inspiree des bornes finales article.
    # Les composantes ±Inf du papier ne sont pas utilisables telles quelles
    # dans Dionysos ; on les remplace ici par une fenetre finie raisonnable
    # basee sur q_abs_bound.
    _T_ = UT.HyperRectangle(
        SVector(
            -0.05,
            -0.05,
            pi / 8,
            params.phi_min,
            0.0,
            -2.0,
            -params.q_abs_bound,
        ),
        SVector(
            0.05,
            0.05,
            pi / 2,
            params.phi_max,
            2.0,
            0.0,
            params.q_abs_bound,
        ),
    )

    _U_ = UT.HyperRectangle(SVector(-params.phi_dot_max), SVector(params.phi_dot_max))
    sys = system(_X_; _U_ = _U_, params = params)

    # Comme pour le benchmark vehicule, on garde une factory tres simple :
    # c'est une instanciation benchmark inspiree de l'article, adaptee a un
    # cadre HyperRectangle fini compatible avec Dionysos.
    return PB.OptimalControlProblem(
        sys,
        _I_,
        _T_,
        state_cost,
        transition_cost,
        terminal_cost,
    )
end

################################################
############ Simple Controllers ################
################################################

function get_constant_controller(u_const)
    u_vec = u_const isa Number ? SVector(u_const) : u_const
    return ST.ConstantController(u_vec)
end

_wrap_to_pi(α) = mod(α + pi, 2pi) - pi

function get_pitch_regulation_controller(θ_ref; kφ = 2.0, kq = 0.5, umax = 13.0)
    f = x -> begin
        θ = x[3]
        φ = x[4]
        q = x[7]

        # Controleur volontairement simple :
        # on fabrique une consigne d'elevator proportionnelle a l'erreur de pitch,
        # puis on commande sa vitesse avec une saturation franche.
        eθ = _wrap_to_pi(θ_ref - θ)
        φref = clamp(kφ * eθ - kq * q, -pi / 3, pi / 8)
        φdot_cmd = clamp(4.0 * (φref - φ), -umax, umax)
        return SVector(φdot_cmd)
    end
    return ST.BlackBoxContinuousController(f)
end

################################################
########## Visualization tools #################
################################################

Base.@kwdef struct DrawParams{T}
    fuselage_length::T
    fuselage_width::T
    wing_chord_draw::T
    wing_span_draw::T
    tail_length::T
    tail_width::T
    elevator_length::T
    elevator_width::T
end

function DrawParams(
    p::Params;
    fuselage_length = max(1.15 * (p.l + p.le), 0.28),
    fuselage_width = 0.035,
    wing_chord_draw = max(0.7 * (p.l + p.le), 0.16),
    wing_span_draw = 0.025,
    tail_length = max(p.l - p.lw, 0.12),
    tail_width = 0.018,
    elevator_length = max(2.0 * p.le, 0.08),
    elevator_width = 0.020,
)
    T = promote_type(
        typeof(fuselage_length),
        typeof(fuselage_width),
        typeof(wing_chord_draw),
        typeof(wing_span_draw),
        typeof(tail_length),
        typeof(tail_width),
        typeof(elevator_length),
        typeof(elevator_width),
    )

    return DrawParams{T}(
        T(fuselage_length),
        T(fuselage_width),
        T(wing_chord_draw),
        T(wing_span_draw),
        T(tail_length),
        T(tail_width),
        T(elevator_length),
        T(elevator_width),
    )
end

rot2(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

# Rectangle centre en c, oriente par θ, de taille (L, W).
function rect_poly(c::SVector{2, T}, θ, L, W) where {T <: Real}
    R = rot2(θ)
    hl, hw = L / 2, W / 2
    pts = (
        c + R * SVector(hl, hw),
        c + R * SVector(hl, -hw),
        c + R * SVector(-hl, -hw),
        c + R * SVector(-hl, hw),
        c + R * SVector(hl, hw),
    )
    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]
    return xs, ys
end

function aircraft_keypoints(p, x, dp::DrawParams)
    px, pz, θ, φ = x[1], x[2], x[3], x[4]

    COM = SVector(Float64(px), Float64(pz))
    e_body = SVector(cos(θ), sin(θ))
    e_elev = SVector(cos(θ + φ), sin(θ + φ))

    nose = COM + 0.5 * dp.fuselage_length * e_body
    tail = COM - 0.5 * dp.fuselage_length * e_body
    wing_center = COM - p.lw * e_body
    elevator_hinge = COM - p.l * e_body
    elevator_center = elevator_hinge - p.le * e_elev

    return COM, nose, tail, wing_center, elevator_hinge, elevator_center, θ, θ + φ
end

function draw_glider!(
    plt,
    p,
    dp::DrawParams,
    x,
    u;
    show_axes = true,
    show_heading = true,
    show_elevator = true,
)
    COM, nose, tail, wing_center, elevator_hinge, elevator_center, θ_body, θ_elev =
        aircraft_keypoints(p, x, dp)

    # Fuselage.
    fuselage_center = 0.5 * (nose + tail)
    fx, fz = rect_poly(fuselage_center, θ_body, dp.fuselage_length, dp.fuselage_width)
    plot!(plt, fx, fz; lw = 1, fill = (true, 0.12), color = :black, label = false)

    # Aile principale, representee simplement comme une petite plaque alignee avec le fuselage.
    wx, wz = rect_poly(wing_center, θ_body, dp.wing_chord_draw, dp.wing_span_draw)
    plot!(plt, wx, wz; lw = 1, fill = (true, 0.10), color = :royalblue, label = false)

    # Segment de queue entre fuselage et charniere.
    plot!(
        plt,
        [tail[1], elevator_hinge[1]],
        [tail[2], elevator_hinge[2]];
        lw = 2,
        color = :black,
        label = false,
    )

    if show_elevator
        ex, ez = rect_poly(elevator_center, θ_elev, dp.elevator_length, dp.elevator_width)
        plot!(plt, ex, ez; lw = 1, fill = (true, 0.18), color = :firebrick, label = false)
        scatter!(plt, [elevator_hinge[1]], [elevator_hinge[2]]; ms = 2, color = :black, label = false)
    end

    if show_axes
        e_body = SVector(cos(θ_body), sin(θ_body))
        n_body = SVector(-sin(θ_body), cos(θ_body))
        axis_scale = 0.25 * dp.fuselage_length

        plot!(
            plt,
            [COM[1], COM[1] + axis_scale * e_body[1]],
            [COM[2], COM[2] + axis_scale * e_body[2]];
            lw = 1.5,
            color = :darkgreen,
            label = false,
        )
        plot!(
            plt,
            [COM[1], COM[1] + axis_scale * n_body[1]],
            [COM[2], COM[2] + axis_scale * n_body[2]];
            lw = 1.2,
            ls = :dot,
            color = :darkorange,
            label = false,
        )
    end

    if show_heading
        plot!(
            plt,
            [COM[1], nose[1]],
            [COM[2], nose[2]];
            lw = 2,
            color = :black,
            label = false,
        )
    end

    return plt
end

function plot_perch!(plt; perch_pos = SVector(0.0, 0.0), radius = 0.05)
    perch = SVector{2, Float64}(perch_pos)
    ts = range(0, 2pi; length = 80)
    xs = [perch[1] + radius * cos(t) for t in ts]
    zs = [perch[2] + radius * sin(t) for t in ts]
    plot!(plt, xs, zs; lw = 2, color = :brown, label = false)
    scatter!(plt, [perch[1]], [perch[2]]; ms = 3, color = :brown, label = false)
    return plt
end

function live_glider_progression(
    p,
    dp,
    x_traj::ST.Trajectory,
    u_traj::ST.Trajectory,
    xl,
    zl;
    domain = nothing,
    every = 1,
    dt = 0.05,
    giffile::Union{Nothing, String} = nothing,
    fps::Int = 20,
    title::Union{Nothing, String} = nothing,
    perch_pos = SVector(0.0, 0.0),
    perch_radius = 0.05,
)
    states = x_traj.seq
    inputs = u_traj.seq

    xs = [x[1] for x in states]
    zs = [x[2] for x in states]
    plot_title = title === nothing ? "" : title

    if giffile !== nothing
        anim = @animate for k in 1:every:length(states)
            plt = plot(;
                aspect_ratio = :equal,
                xlims = xl,
                ylims = zl,
                legend = false,
                size = (800, 500),
                title = plot_title,
                xlabel = "x (m)",
                ylabel = "z (m)",
            )
            if domain !== nothing
                plot!(plt, domain; color = :grey, opacity = 0.10)
            end
            plot_perch!(plt; perch_pos = perch_pos, radius = perch_radius)
            plot!(plt, xs, zs; lw = 1.5, color = :steelblue, label = false)
            uk = isempty(inputs) ? SVector(0.0) : ((k <= length(inputs)) ? inputs[k] : inputs[end])
            draw_glider!(plt, p, dp, states[k], uk)
        end

        gif(anim, giffile; fps = fps)
        return anim
    end

    for k in 1:every:length(states)
        plt = plot(;
            aspect_ratio = :equal,
            xlims = xl,
            ylims = zl,
            legend = false,
            size = (800, 500),
            title = plot_title,
            xlabel = "x (m)",
            ylabel = "z (m)",
        )
        if domain !== nothing
            plot!(plt, domain; color = :grey, opacity = 0.10)
        end
        plot_perch!(plt; perch_pos = perch_pos, radius = perch_radius)
        plot!(plt, xs, zs; lw = 1.5, color = :steelblue, label = false)
        uk = isempty(inputs) ? SVector(0.0) : ((k <= length(inputs)) ? inputs[k] : inputs[end])
        draw_glider!(plt, p, dp, states[k], uk)
        display(plt)
        sleep(dt)
    end

    return nothing
end

# Petit alias de compatibilite :
# certains helpers historiques du depot attendent un point d'entree
# `live_vehicle_progression` fourni par les benchmarks vehicule.
# On le redirige simplement vers la version glider.
function live_vehicle_progression(args...; kwargs...)
    return live_glider_progression(args...; kwargs...)
end

end # module
