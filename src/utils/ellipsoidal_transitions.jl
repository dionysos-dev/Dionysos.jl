import MathematicalSystems
MS = MathematicalSystems

import LinearAlgebra as LA

import HybridSystems

AffineSys = Union{
    HybridSystems.NoisyConstrainedAffineControlDiscreteSystem,
    HybridSystems.ConstrainedAffineControlDiscreteSystem,
    HybridSystems.HybridSystems.ConstrainedAffineControlMap,
}

function format_input_set(rec::HyperRectangle)
    n = get_dims(rec)
    Uaux = LA.diagm(1:n)
    U = [(Uaux .== i) ./ rec.ub[i] for i in 1:n]
    return U
end

function format_input_set(elli::Ellipsoid)
    return [get_root(elli)]
end

function format_input_set(iset::LazySetIntersection)
    result = Any[]
    for set in iset.sets
        Base.append!(result, format_input_set(set))
    end
    return result
end

function format_noise_set(rec::HyperRectangle)
    return get_vertices(rec)
end

function get_controller_matrices(m)
    dims = size(m)
    nu = dims[1]
    return m[:, 1:(end - 1)], SVector{nu}(m[:, end])
end

function _getμν(L, subsys::AffineSys)
    n_x = length(subsys.c)
    return (
        LazySets.vertices_list(IA.IntervalBox(IA.interval(-x, x) for x in L[1:n_x])),
        LazySets.vertices_list(IA.IntervalBox(subsys.D * subsys.W...)),
    )
end

function hasTransition(
    c,
    u,
    Ep::Ellipsoid,
    subsys::AffineSys,
    L,
    S,
    U,
    maxRadius,
    maxΔu,
    optimizer;
    λ = 0.01,
)
    Pp = Ep.P
    cp = Ep.c

    eye(n) = LA.diagm(ones(n))
    A = subsys.A
    B = subsys.B
    g = subsys.c
    n = length(g)
    m = size(U[1], 2)

    μ, ν = _getμν(L, subsys)
    N_μ = length(μ)
    N_ν = length(ν)
    N_u = length(U)

    model = Model(optimizer)
    @variable(model, C[i = 1:n, j = 1:n])
    @variable(model, X[i = 1:n, j = 1:n], PSD)
    @variable(model, F[i = 1:m, j = 1:n])
    @variable(model, ell[i = 1:m, j = 1:1])
    @variable(model, bta[i = 1:N_μ, j = 1:N_ν] >= 0)
    @variable(model, tau[i = 1:N_u] >= 0)
    @variable(model, gamma >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, r >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϵ >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A * C + B * F
        gt, g + B * ell
    end)

    t(x) = transpose(x)

    z = zeros(n, 1)

    for i in 1:N_μ
        for j in 1:N_ν
            aux = @expression(
                model,
                A * hcat(c) + hcat(gt) - hcat(cp) +
                hcat(Vector(μ[i])) * (r + δu) +
                hcat(Vector(ν[j]))
            )#
            @constraint(
                model,
                [
                    bta[i, j]*eye(n) z t(At)
                    t(z) 1-bta[i, j] t(aux)
                    At aux inv(Pp)
                ] >= eye(2 * n + 1) * 1e-4,
                PSDCone()
            )
        end
    end

    for i in 1:N_u
        n_ui = size(U[i], 1)
        @constraint(
            model,
            [
                tau[i]*eye(n) z t(U[i] * F)
                t(z) 1-tau[i] t(U[i] * ell)
                U[i]*F U[i]*ell eye(n_ui)
            ] >= eye(n + n_ui + 1) * 1e-4,
            PSDCone()
        )
    end
    n_S = size(S, 1)
    @constraint(
        model,
        [
            gamma*eye(n) z [t(C) t(F) z]*t(S)
            t(z) J-gamma [t(c) t(ell) 1]*t(S)
            S*t([t(C) t(F) z]) S*t([t(c) t(ell) 1]) eye(n_S)
        ] >= eye(n + n_S + 1) * 1e-4,
        PSDCone()
    )

    u = hcat(u)
    @constraint(
        model,
        [
            ϕ*eye(n) z t(F)
            t(z) δu-ϕ t(ell - u)
            (F) (ell-u) eye(m)
        ] >= eye(n + m + 1) * 1e-4,
        PSDCone()
    )

    @constraint(model, [
        eye(n) t(C)
        C r*eye(n)
    ] >= eye(n * 2) * 1e-4, PSDCone())

    # @constraint(model,diag(C).>=ones(n,1)*0.01)
    @constraint(model, r <= maxRadius^2)
    @constraint(model, δu <= maxΔu^2)
    @constraint(model, diag(C) .>= ones(n, 1) * ϵ)

    @objective(model, Min, -ϵ + λ * J)# -tr(C)) #TODO regularization ? 

    optimize!(model)
    if solution_summary(model).termination_status == MOI.OPTIMAL
        C = value.(C)
        El = Ellipsoid(transpose(C) \ eye(n) / C, c)
        kappa = [value.(F) / (C) value.(ell)]
        cost = value(J)
    else
        El = nothing
        kappa = nothing
        cost = nothing
    end
    # println("$(solution_summary(model).solve_time) s")
    return El, kappa, cost
end

#     _has_transition(A, B, g, U, W, L, c, P, cp, Pp, optimizer)

# Verifies whether a controller u(x)=K(x-c)+ell exists for `subsys` satisfying input requirements
# defined by `U` ans performs a sucessful transitions from a starting set Bs = {(x-c)'P(x-c) ≤ 1}
# to the final set Bs = {(x-c)'P(x-c) ≤ 1}. A tight upper bound on the transition cost c(x,u) is 
# minimized where c(x,u) = |L*[x; u; 1]|^2, from the parameter `L` and each columm of the matrix  
# `W` defines a vertex of the polytope from which additive disturbance are drawn. 

# The input restrictions are defined by the list U as |U[i]*u| ≤ 1. `optimizer` must be a JuMP 
# SDP optimizer (e.g., Mosek, SDPA, COSMO, ...).

# ## Note
# This implements the optimization problem presented in Corollary 1 of the following paper 
# https://arxiv.org/pdf/2204.00315.pdf
function _has_transition(A, B, g, U, W, L, c, P, cp, Pp, optimizer)
    eye(n) = LA.diagm(ones(n))
    n = length(c)
    m = size(U[1], 2)
    N = size(W, 2)
    p = size(W, 1)
    Nu = length(U)

    model = Model(optimizer)
    @variable(model, K[i = 1:m, j = 1:n])
    @variable(model, ell[i = 1:m, j = 1:1])
    @variable(model, bta[i = 1:N] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, gamma >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A + B * K
        gt, g + B * ell
    end)

    t(x) = transpose(x)

    z = zeros(n, 1)

    for i in 1:N
        w = W[:, i]
        aux = A * hcat(c) + hcat(gt) - hcat(cp) + hcat(w)
        @constraint(
            model,
            [
                bta[i]*P z t(At)
                t(z) 1-bta[i] t(aux)
                At aux inv(Pp)
            ] >= eye(2 * n + 1) * 1e-4,
            PSDCone()
        )
    end

    for i in 1:Nu
        n_ui = size(U[i], 1)
        @constraint(
            model,
            [
                tau[i]*P z t(U[i] * K)
                t(z) 1-tau[i] t(U[i] * ell)
                U[i]*K U[i]*ell eye(n_ui)
            ] >= eye(n + n_ui + 1) * 1e-4,
            PSDCone()
        )
    end
    """n_S = size(L, 1)
    @constraint(
        model,
        [
            gamma*P z [LA.I t(K) z]*t(L)
            t(z) J-gamma [t(c) t(ell) 1]*t(L)
            t([LA.I t(K) z] * t(L)) t([t(c) t(ell) 1] * t(L)) eye(n_S)
        ] >= eye(n + n_S + 1) * 1e-4,
        PSDCone()
    )"""
    n  = length(c)
    m  = size(U[1], 2)
    nS = size(L, 1)

    I_n = Matrix{Float64}(LA.I, n, n)
    Z_n1 = zeros(n, 1)

    # This is n × (n+m+1)
    G = hcat(I_n, transpose(K), Z_n1)

    # This is 1 × (n+m+1)
    d = hcat(transpose(c), transpose(ell), 1.0)
    @assert size(L, 2) == n + m + 1 "Expected L to have $(n+m+1) columns (for [x;u;1]), got size(L)=$(size(L))."
    @constraint(model,
        [
            gamma * P                zeros(n,1)          G * transpose(L)
            zeros(1,n)               J - gamma           d * transpose(L)
            L * transpose(G)         L * transpose(d)    Matrix{Float64}(LA.I, nS, nS)
        ] in PSDCone()
    )

    @objective(model, Min, J)

    optimize!(model)

    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    kappa = [value.(K) value.(ell)]
    cost = value(J)
    # println("$(solution_summary(model).solve_time) s")
    return ans, kappa, cost
end

function _has_transition(
    affsys::AffineSys,
    E1::Ellipsoid,
    E2::Ellipsoid,
    U,
    W,
    S,
    optimizer,
)
    ans, kappa, cost = _has_transition(
        affsys.A,
        affsys.B,
        affsys.c,
        U,
        W,
        S,
        E1.c,
        E1.P,
        E2.c,
        E2.P,
        optimizer,
    )
    K, ℓ = get_controller_matrices(kappa)
    cont = MS.AffineMap(K, ℓ-K*E1.c)
    return ans, cont, cost
end

#
#    _compute_base_cell(r::SVector{S})
#
# Computes a polyhedron containing the base hyperrectangular cell, centered at the origin
# and with the i-th side lenght given by `2*r[i]`. 
#
function _compute_base_cell(r::SVector{S}) where {S}
    baseCellList = []
    for i in 1:S
        vec = SVector{S}(1:S .== i)
        append!(
            baseCellList,
            [Polyhedra.HalfSpace(-vec, r[i]) ∩ Polyhedra.HalfSpace(vec, r[i])],
        )
    end
    return Polyhedra.polyhedron(intersect(baseCellList...))
end

#
#    _provide_P(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem, optimizer)
#
# If `subsys` is a stabilizable system, finds the matrix `P` and the state-feedback gain `K`
# that satisfy the discrete-time Lyapunov inequality (A+BK)'P(A+BK)-P < 0. The condition number
# of `P` is minimized. `optimizer` must be a JuMP SDP optimizer.

function _provide_P(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem, optimizer)
    eye(n) = LA.diagm(ones(n))
    A = subsys.A
    B = subsys.B
    n = size(A, 1)
    m = size(B, 2)

    model = Model(optimizer)
    @variable(model, L[i = 1:m, j = 1:n])
    @variable(model, S[i = 1:n, j = 1:n], PSD)
    @variable(model, gamma)

    t(x) = transpose(x)

    @constraint(model, [
        S t(A * S + B * L)
        A * S+B * L S
    ] >= 1e-4 * eye(2n), PSDCone())
    @constraint(model, eye(n) >= S, PSDCone())
    @constraint(model, S >= gamma * eye(n), PSDCone())

    @objective(model, Max, gamma)

    optimize!(model)

    P = inv(value.(S))
    K = value.(L) * P
    gamma = value(gamma)
    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    return ans, K, P, gamma
end

###############################################################################################################################

# U: each entry of U is a quadratic constraint on the input
# W: each column is vertex of the polytopic noise 
# x(k+1) = Ax(k)+Bu(k)+g+Ew(k), with u(k)∈E(0,U^T U), w(k)∈W
# W[:,i] = vertex i of the polytop
function transition_fixed(A, B, c, D, U, W, S, c1, P1, c2, P2, optimizer)
    W = D * W
    eye(n) = LA.diagm(ones(n))
    nx = length(c) # dimension of the state
    nu = size(U[1], 2) # dimension of the input
    nw = size(W, 1) # dimension of the noise (not use for now, we could add matrix of the noise)
    Nu = length(U) # number of quadratic input constraints
    Nw = size(W, 2) # number of vertex of the polytopic noise

    model = Model(optimizer)
    @variable(model, K[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A + B * K
        ct, c + B * ell
    end)

    t(x) = transpose(x)
    z = zeros(nx, 1)

    for i in 1:Nw
        w = W[:, i]
        aux = A * hcat(c1) + hcat(ct) - hcat(c2) + hcat(w)
        @constraint(
            model,
            [
                beta[i]*P1 z t(At)
                t(z) 1-beta[i] t(aux)
                At aux inv(P2)
            ] >= eye(2 * nx + 1) * 1e-4,
            PSDCone()
        )
    end

    for i in 1:Nu
        n_ui = size(U[i], 1)
        @constraint(
            model,
            [
                tau[i]*P1 z t(U[i] * K)
                t(z) 1-tau[i] t(U[i] * ell)
                U[i]*K U[i]*ell eye(n_ui)
            ] >= eye(nx + n_ui + 1) * 1e-4,
            PSDCone()
        )
    end
    n_S = size(S, 1)
    @constraint(
        model,
        [
            γ*P1 z [LA.I t(K) z]*t(S)
            t(z) J-γ [t(c1) t(ell) 1]*t(S)
            S*t([LA.I t(K) z]) S*t([t(c1) t(ell) 1]) eye(n_S)
        ] >= eye(nx + n_S + 1) * 1e-4,
        PSDCone()
    )

    @objective(model, Min, J)

    optimize!(model)

    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    kappa = [value.(K) value.(ell)]
    cost = value(J)

    return ans, kappa, cost
end

function transition_fixed(
    affsys::AffineSys,
    E1::Ellipsoid,
    E2::Ellipsoid,
    U,
    W,
    S,
    optimizer,
)
    ans, kappa, cost = transition_fixed(
        affsys.A,
        affsys.B,
        affsys.c,
        affsys.D,
        U,
        W,
        S,
        E1.c,
        E1.P,
        E2.c,
        E2.P,
        optimizer,
    )
    K, ℓ = get_controller_matrices(kappa)
    cont = MS.AffineMap(K, ℓ-K*E1.c)
    return ans, cont, cost
end

function _getμν(L, nx, D, W)
    vertices_matrix = D * W
    noise_vertices = [vertices_matrix[:, i] for i in 1:size(vertices_matrix, 2)]
    return (
        LazySets.vertices_list(IA.IntervalBox(IA.interval(-x, x) for x in L[1:nx])),
        noise_vertices,
    )
end

# the dynamic: Ax+Bu+c+Dw
# linearization point: (̄x,̄u, w) = (c,u,0)
# inputs constraint: u
# polytopic noise: W (ach column is a vertex)
# objective function: S
function transition_backward(
    A,
    B,
    c,
    D,
    c2,
    P2,
    c1,
    u,
    U,
    W,
    S,
    Lip,
    optimizer;
    maxδx = maxδx,
    maxδu = maxδu,
    λ = 0.01,
)
    eye(n) = LA.diagm(ones(n))
    nx = length(c) #dimension of the state
    nu = size(U[1], 2) #dimension of the input
    μ, ν = _getμν(Lip, nx, D, W)
    Nx = length(μ) #number of vertex of the hyperrectangle: 2^nx
    Nw = length(ν) #number of vertex of the polytopic noise
    Nu = length(U) #number of constraints on u

    model = Model(optimizer)
    @variable(model, L[i = 1:nx, j = 1:nx], PSD)
    @variable(model, F[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx, j = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, δx >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A * L + B * F
        ct, c + B * ell
    end)

    t(x) = transpose(x)
    z = zeros(nx, 1)

    for i in 1:Nx
        for j in 1:Nw
            aux = @expression(
                model,
                A * hcat(c1) + hcat(ct) - hcat(c2) +
                hcat(Vector(μ[i])) * (δx + δu) +
                hcat(Vector(ν[j]))
            )
            @constraint(
                model,
                [
                    beta[i, j]*eye(nx) z t(At)
                    t(z) 1-beta[i, j] t(aux)
                    At aux inv(P2)
                ] >= eye(2 * nx + 1) * 1e-4,
                PSDCone()
            )
        end
    end

    for i in 1:Nu
        n_ui = size(U[i], 1)
        @constraint(
            model,
            [
                tau[i]*eye(nx) z t(U[i] * F)
                t(z) 1-tau[i] t(U[i] * ell)
                U[i]*F U[i]*ell eye(n_ui)
            ] >= eye(nx + n_ui + 1) * 1e-4,
            PSDCone()
        )
    end
    n_S = size(S, 1)
    @constraint(
        model,
        [
            γ*eye(nx) z [t(L) t(F) z]*t(S)
            t(z) J-γ [t(c1) t(ell) 1]*t(S)
            S*t([t(L) t(F) z]) S*t([t(c1) t(ell) 1]) eye(n_S)
        ] >= eye(nx + n_S + 1) * 1e-4,
        PSDCone()
    )

    u = hcat(u)
    @constraint(
        model,
        [
            ϕ*eye(nx) z t(F)
            t(z) δu-ϕ t(ell - u)
            (F) (ell-u) eye(nu)
        ] >= eye(nx + nu + 1) * 1e-4,
        PSDCone()
    )

    @constraint(model, [
        eye(nx) t(L)
        L δx*eye(nx)
    ] >= eye(nx * 2) * 1e-4, PSDCone())

    @constraint(model, δx <= maxδx^2)
    @constraint(model, δu <= maxδu^2)

    @variable(model, t)
    u_q = [L[i, j] for j in 1:nx for i in 1:j]

    @constraint(model, vcat(t, 1, u_q) in MOI.LogDetConeTriangle(nx))

    @objective(model, Min, λ * J + (1 - λ) * (-t))
    optimize!(model)

    if solution_summary(model).termination_status == MOI.OPTIMAL
        L = value.(L)
        P = transpose(L) \ eye(nx) / L
        kappa = [value.(F) / (L) value.(ell)]
        cost = value(J)
    else
        P = nothing
        kappa = nothing
        cost = nothing
    end
    return P, kappa, cost
end

function transition_backward(
    affsys::AffineSys,
    E2::Ellipsoid,
    c1,
    u,
    U,
    W,
    S,
    Lip,
    optimizer;
    maxδx = 100,
    maxδu = 10.0 * 2,
    λ = 0.01,
)
    P1, kappa, cost = transition_backward(
        affsys.A,
        affsys.B,
        affsys.c,
        affsys.D,
        E2.c,
        E2.P,
        c1,
        u,
        U,
        W,
        S,
        Lip,
        optimizer;
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
    )
    if P1 !== nothing
        K, ℓ = get_controller_matrices(kappa)
        cont = MS.AffineMap(K, ℓ-K*c1)
        return Ellipsoid(P1, c1), cont, cost
    else
        return nothing, nothing, nothing
    end
end