export build_symbolic_context, backward_step!, run_backward_chain!

import Dionysos
import IntervalArithmetic as IA
import LinearAlgebra as LA

const DI = Dionysos
const ST = DI.System
const SY = DI.Symbolic
const UT = DI.Utils

struct EllipsoidalBackwardContext{P, C, CFG, SYM, TX, TU, TS, TB} # améliorer les typages
    problem::P
    candidate::C
    config::CFG
    symbolic::SYM
    xs::TX
    us::TU
    K::Int
    fsymbolic
    x
    u
    w
    ΔX
    ΔU
    ΔW
    Uformat
    Wformat
    S::TS
    backend::TB
    maxδx::Float64
    maxδu::Float64
    λ::Float64
    terminal_radius::Float64
end

function _identity_transition_cost(nx::Int, nu::Int) # l'utilisateur devrait lui meme fournir la transi cost fonction (pour le moment c'est pratiqque)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

function build_symbolic_context(problem, candidate, config, symbolic_builder) # je pourrais cancel ça et tout condenser
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    K = length(us)
    @assert length(xs) == K + 1

    sym = symbolic_builder(problem, candidate, config)

    fsymbolic = sym.fsymbolic
    x = sym.x
    u = sym.u
    w = sym.w
    ΔX = sym.ΔX
    ΔU = sym.ΔU
    ΔW = sym.ΔW
    Uformat = sym.Uformat
    Wformat = sym.Wformat

    nx = length(xs[1])
    nu = length(us[1])
    S = _identity_transition_cost(nx, nu)

    opts = config.options
    maxδx = Float64(opts.maxδx)
    maxδu = Float64(opts.maxδu)
    λ = Float64(opts.λ)
    terminal_radius = Float64(opts.rayon_terminal)

    return EllipsoidalBackwardContext(
        problem,
        candidate,
        config,
        sym,
        xs,
        us,
        K,
        fsymbolic,
        x,
        u,
        w,
        ΔX,
        ΔU,
        ΔW,
        Uformat,
        Wformat,
        S,
        config.backend,
        maxδx,
        maxδu,
        λ,
        terminal_radius,
    )
end

function backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    xk = collect(ctx.xs[k])
    uk = collect(ctx.us[k])
    wk = collect(zeros(length(ctx.w)))

    Xbar = IA.IntervalBox(xk .+ ctx.ΔX)
    Ubar = IA.IntervalBox(uk .+ ctx.ΔU)
    Wbar = IA.IntervalBox(wk .+ ctx.ΔW)

    affineSys, L = ST.buildAffineApproximation(
        ctx.fsymbolic,
        ctx.x,
        ctx.u,
        ctx.w,
        xk,
        uk,
        wk,
        Xbar,
        Ubar,
        Wbar,
    )
    # println("\n\nk=", k, " | ū=", uk, " | Lips=", L)


    E_prev, kappa, cost = UT.transition_backward(
        affineSys,
        E_next,
        xk,
        uk,
        ctx.Uformat,
        ctx.Wformat,
        ctx.S,
        L,
        ctx.backend;
        maxδx = ctx.maxδx,
        maxδu = ctx.maxδu,
        λ = ctx.λ,
    )

    if E_prev === nothing || kappa === nothing
        println("My transi is impossible")
        return BackwardStepRecord{Any, Any, NamedTuple}(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (; L),
        )
    end

    return BackwardStepRecord{Any, Any, NamedTuple}(
        k,
        :ok,
        Float64(cost),
        E_prev,
        kappa,
        (; L),
    )
end

function run_backward_chain!(ctx::EllipsoidalBackwardContext)
    t0 = time()

    nx = length(ctx.xs[end])
    PN = Matrix{Float64}(LA.I, nx, nx) * (1.0 / ctx.terminal_radius^2) # je sais pas si je vais changer ou non cette logique
    E_next = UT.Ellipsoid(PN, collect(ctx.xs[end]))
    

    steps = BackwardStepRecord{Any, Any, NamedTuple}[]
    ellipsoids = UT.Ellipsoid[E_next]
    kappas = Any[]

    for k in ctx.K:-1:1
        rec = backward_step!(ctx, k, E_next)
        push!(steps, rec)
       
        if rec.status == :infeasible
            return EllipsoidalCertificationResult(
                false,
                k,
                Float64(time() - t0),
                steps,
                nothing,
                (; ellipsoids, kappas),
            )
        end
        
        #println("[",k,"] is a succes\n")
        E_next = rec.ellipsoid
        
        #println("the volume is ",UT.get_volume(E_next),"\n\n")
        push!(ellipsoids, rec.ellipsoid)
        push!(kappas, rec.kappa)
        #println(E_next,"\n\n")
    end

    return EllipsoidalCertificationResult(
        true,
        nothing,
        Float64(time() - t0),
        steps,
        kappas,
        (; ellipsoids, kappas),
    )
end
