# The half of `Dionysos.Utils.PathCompleteFramework` that computes certificates, kept here for the
# research scripts. Dionysos now leaves that job to PathCompleteCertificates.jl, whose certificates
# a package extension turns into `PathCompleteFramework.PCLF`. The functions below are the former
# library code, verbatim: the experiments of `research/HSCC2027` and `research/BisimulationQuotient`
# keep calling exactly the programs their published numbers came from. A different formulation of
# the same template finds another certificate among the many optimal ones, and with it other level
# sets, other cells and other counts.
module PCLFProducer

import Dionysos
const UT = Dionysos.Utils
import Dionysos.Utils.PathCompleteFramework:
    LabDigraph, edgeList_to_LabDigraph, AbstractPiece, PolyhedralPiece, PCLF
import HybridSystems
import JuMP
const MOI = JuMP.MOI
import LinearAlgebra

function generate_DeBruijn_edges(M::Int, k::Int; dual::Bool = false)
    @assert M ≥ 1
    @assert k ≥ 0

    # Common Lyapunov function graph:
    if k == 0
        edges = Vector{Tuple{Int, Int, Int}}()
        v = (1)
        for s in 1:M
            push!(edges, (v, v, s))
        end
        return edgeList_to_LabDigraph(edges)
    end

    # Otherwise: 
    edges = Vector{Tuple{NTuple{k, Int}, NTuple{k, Int}, Int}}()

    iterables = ntuple(_ -> 1:M, k)
    nodes = collect(Iterators.product(iterables...))

    for node in nodes
        for state in nodes
            if node[2:end] == state[1:(end - 1)]
                if !dual
                    push!(edges, (node, state, state[end]))
                else
                    push!(edges, (state, node, state[end]))
                end
            end
        end
    end

    return edgeList_to_LabDigraph(edges)
end

"""
Compute a path-complete Lyapunov function (PCLF) with **symmetric polyhedral pieces
having 2n faces** for a switched linear system.

Each node `s` of the graph is associated with a polyhedral Lyapunov function of the form:

    V_s(x) = max_i |(G_s x)_i| / w_s[i]

whose sublevel sets are polytopes:

    { x : -γ w_s ≤ G_s x ≤ γ w_s }.

# Method
The method constructs and solves a feasibility linear program (LP) using bisection on γ.
For each edge (u → v, σ), it enforces:

    |G_v A_σ G_u^{-1}| * w_u ≤ γ w_v,

where the absolute value is taken elementwise.

# Arguments
- `f`: hybrid system containing the system matrices `A_σ`
- `D`: labeled directed graph defining the PCLF structure
- `optimizer`: JuMP optimizer

# Keyword arguments
- `Gmats`: choice of matrices G_s (identity, Dict, or Vector)
- `tol`: tolerance for bisection on γ
- `maxiter`: maximum number of bisection iterations
- `MLF`: if true, extracts the Lyapunov pieces
- `verbose`: enable solver output
- `min_w`: lower bound to enforce strict positivity of w

# Returns
- `PCLF`: structure containing the graph, Lyapunov pieces, and JSR approximation

# Notes
- The resulting Lyapunov functions are structured and correspond to
  weighted ∞-norms in transformed coordinates.
- This approach is computationally efficient but may be conservative.

# Reference
This symmetric `2n`-face construction of polyhedral path-complete Lyapunov
functions follows [athanasopoulos2019polyhedral](@cite).
"""
function compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f::HybridSystems.HybridSystem,
    D::LabDigraph,
    optimizer;
    Gmats = :identity,
    tol = 1e-5,
    maxiter = 100,
    MLF = false,
    verbose = false,
    min_w = 1e-3,
)

    # --- extract matrices from resetmaps ---
    A = UT.mode_matrices(f)

    # --- vertices and indexing for w variables ---
    verts = collect(D.verts)
    l_s = length(verts)
    index_of = Dict{typeof(verts[1]), Int}()
    for (i, v) in enumerate(verts)
        index_of[v] = i
    end

    # Normalize Gmats input to a Vector indexed 1..l_s
    # --- Build G_by_idx ---
    n = size(A[1], 1)
    G_by_idx = Vector{Matrix{Float64}}(undef, l_s)

    if Gmats === :identity
        # Default: all G_s = I
        for i in 1:l_s
            G_by_idx[i] = Matrix{Float64}(LinearAlgebra.I, n, n)
        end

    elseif isa(Gmats, Dict)
        for (i, v) in enumerate(verts)
            @assert haskey(Gmats, v) "Gmats missing vertex $v"
            G_by_idx[i] = Array(Gmats[v])
        end

    elseif isa(Gmats, AbstractVector)
        @assert length(Gmats) >= l_s "Gmats vector too short"
        for i in 1:l_s
            G_by_idx[i] = Array(Gmats[i])
        end

    else
        error("Gmats must be :identity, Dict or Vector of matrices")
    end

    # sizes and inverses
    n = size(G_by_idx[1], 1)
    for i in 1:l_s
        @assert size(G_by_idx[i], 1) == n && size(G_by_idx[i], 2) == n "All G matrices must be n×n"
    end
    Ginv = [inv(G_by_idx[i]) for i in 1:l_s]

    # --- Precompute M matrices for each edge: M = |G_d * A_sigma * G_s^{-1}| ---
    edge_list = D.edges  # expected iterable of (u,v,label)
    Mlist = Vector{Tuple{Int, Int, Int, Matrix{Float64}}}()  # (ui,vi,sigma,M)
    for (u, v, label) in edge_list
        ui = index_of[u]
        vi = index_of[v]
        σ = Int(label)
        M = abs.(G_by_idx[vi] * A[σ] * Ginv[ui])   # nonnegative n×n matrix
        push!(Mlist, (ui, vi, σ, M))
    end

    # --- initial upper bound b: max row sum among the M matrices ---
    #
    # It must be taken on the `M = |G_v A_σ G_u^{-1}|`, NOT on the raw `A_σ`. The constraint the
    # bisection tests is `M w_u <= γ w_v`, so giving every `w` the same value makes it feasible
    # exactly at `γ = max row sum of M` — a bracket the search is *guaranteed* to start inside.
    # `opnorm(A_σ, Inf)` is that same quantity only when the templates are the identity; under a
    # rotated template it can be strictly smaller, and then the top of the bracket is itself
    # infeasible and a certifiable system is reported as having no certificate.
    a = 0.0
    b = 0.0
    for (_, _, _, M) in Mlist
        b = max(b, maximum(sum(M; dims = 2)))
    end

    # See `compute_quadratic_pieces_pclf` for why the lower bracket is language-aware.
    for (u, v, label) in D.edges
        u == v || continue
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(A[Int(label)]))))
    end

    # Whether one trial rate admits piece weights, as a closure so the bracket can be tested at `b`
    # itself rather than only at points strictly inside it.
    function feasible_at_rate(gamma)
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        # variables: w_i in R^n with strict positivity lower bound min_w
        wvars = [JuMP.@variable(model, [1:n]) for i in 1:l_s]
        for i in 1:l_s, k in 1:n
            JuMP.@constraint(model, wvars[i][k] >= min_w)
        end

        # constraints: for every precomputed M (ui->vi): M * w_ui <= gamma * w_vi
        for (ui, vi, _, M) in Mlist
            for k in 1:n
                JuMP.@constraint(
                    model,
                    sum(M[k, p] * wvars[ui][p] for p in 1:n) <= gamma * wvars[vi][k]
                )
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        return st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT
    end

    # --- bisection ---
    #
    # The top of the bracket is TESTED, not assumed feasible. Without this the search reports failure
    # whenever the optimum sits exactly at `b`, because the lower bracket `a = ρ(A_m)` over self-loops
    # then equals `b`, `b - a > tol` is false at once, the loop body never runs and `any_feasible`
    # stays false. That is not a corner case: `A = diag(0.6, 0.48)` under a π/6-rotated template has
    # Perron root exactly 0.6 — a similarity transform preserves the spectrum — so `a` and `b` meet
    # on the true answer and a perfectly good certificate was reported as absent.
    any_feasible = feasible_at_rate(b)

    iter = 0
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        gamma = (a + b) / 2
        if feasible_at_rate(gamma)
            any_feasible = true
            b = gamma
        else
            a = gamma
        end
    end

    # With a fixed template and free positive weights a feasible rate ALWAYS exists — give every `w`
    # the same value and the constraints hold at `γ = max row sum of M`, which is how `b` is chosen.
    # So reaching here with nothing feasible does not mean "no certificate"; it means the solver could
    # not certify even the bracket that is feasible by construction, which is a numerical failure and
    # is reported as such. `Inf` is kept for that case rather than returning a rate no solve backed.
    if !any_feasible
        @warn "compute_symmetric_2n_faces_polyhedral_pieces_pclf: the solver rejected even the \
               guaranteed bracket γ = $(b), which is feasible by construction (equal weights). This \
               is a numerical failure, not an absent certificate. Returning JSRapprox = Inf."
        return PCLF(D, Dict{Any, AbstractPiece}(), Inf)
    end

    gamma = b

    # A rate is always returned, and it is the caller's job to ask whether it contracts. Saying so
    # here because the previous version conflated "no contraction" with "no answer", and a rate of,
    # say, 1.2 is strictly more useful than `Inf`: it says how far the template is from working.
    if gamma >= 1
        @warn "compute_symmetric_2n_faces_polyhedral_pieces_pclf: best rate for this graph and \
               template is $(gamma) >= 1, so it certifies no contraction." maxlog = 1
    end

    # --- final solve to extract pieces (if requested) ---
    pieces = Dict{Any, AbstractPiece}()
    if MLF
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        wvars = [JuMP.@variable(model, [1:n]) for i in 1:l_s]

        # enforce strict positivity
        for i in 1:l_s, k in 1:n
            JuMP.@constraint(model, wvars[i][k] >= min_w)
        end

        for (ui, vi, σ, M) in Mlist
            for k in 1:n
                JuMP.@constraint(
                    model,
                    sum(M[k, p] * wvars[ui][p] for p in 1:n) <= gamma * wvars[vi][k]
                )
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        if !(st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT)
            @warn "Final polyhedral LP not feasible/optimal; status = $st"
        else
            for (i, v) in enumerate(verts)
                wnum = JuMP.value.(wvars[i])
                # numeric safety: enforce positivity
                wnum .= max.(wnum, min_w)
                pieces[v] = PolyhedralPiece(G_by_idx[i], Array(wnum))
            end
        end
    end

    return PCLF(D, pieces, gamma)
end

"""
Compute a path-complete Lyapunov function (PCLF) with **general polyhedral pieces**
defined over a partition of the state space into cones.

Each node `s` is associated with a piecewise-linear Lyapunov function:

    V_s(x) = max_i |p_{s,i}ᵀ x|,

where the rows of a matrix `P_s` define the supporting hyperplanes of the polytope.

# Method
The method formulates a feasibility linear program (LP) based on:

1. Positivity constraints ensuring V_s(x) ≥ 0 on each cone
2. Dominance constraints ensuring correct piecewise structure
3. Decrease conditions along edges:

       V_v(A_σ x) ≤ ρ V_u(x)

These constraints are enforced on the extreme rays of the cones in `partitions`.

A bisection on ρ is used to approximate the joint spectral radius (JSR).

# Arguments
- `f`: hybrid system containing the system matrices `A_σ`
- `D`: labeled directed graph defining the PCLF structure
- `optimizer`: JuMP optimizer
- `partitions`: dictionary mapping each node to a list of cones (matrices of rays)

# Keyword arguments
- `tol`: tolerance for bisection on ρ
- `maxiter`: maximum number of iterations
- `MLF`: if true, extracts the Lyapunov pieces
- `verbose`: enable solver output
- `min_c`: lower bound on auxiliary scalar variables

# Returns
- `PCLF`: structure containing the graph, Lyapunov pieces, and JSR approximation

# Notes
- This method allows for more general polyhedral Lyapunov functions than the
  symmetric 2n-face construction.
- The number of faces depends on the number of rows of `P_s`.
- Less conservative but computationally more expensive.
- The quality depends on the chosen cone partition.
"""
function compute_polyhedral_pieces_pclf(
    f::HybridSystems.HybridSystem,
    D::LabDigraph,
    optimizer,
    partitions;
    tol = 1e-5,
    maxiter = 100,
    MLF = false,
    verbose = false,
    min_c = 1e-3,
)

    # --- extract matrices from resetmaps ---
    A = UT.mode_matrices(f)

    # --- vertices and indexing ---
    verts = collect(D.verts)
    l_s = length(verts)
    index_of = Dict{Any, Int}()
    for (i, v) in enumerate(verts)
        index_of[v] = i
    end

    # --- check partitions and infer dimension ---
    @assert haskey(partitions, verts[1]) "Partitions missing for vertex $(verts[1])"
    @assert !isempty(partitions[verts[1]]) "Node $(verts[1]) must have at least one cone"

    n = size(partitions[verts[1]][1], 1)

    for v in verts
        @assert haskey(partitions, v) "Partitions missing for vertex $v"
        @assert !isempty(partitions[v]) "Node $v must have at least one cone"
        for Xi in partitions[v]
            @assert size(Xi, 1) == n "All cones must live in R^n"
        end
    end

    # --- linear form helper ---
    linrow(P, i, x) = sum(P[i, k] * x[k] for k in 1:n)

    # --- solve feasibility LP for a fixed rho ---
    function feasibility_at_rho(rho::Float64; extract_solution::Bool = false)
        model = JuMP.Model(optimizer)
        if !verbose
            JuMP.set_silent(model)
        end

        # variables: one matrix P_s per node, one c_s per node
        Pvars = Dict{Any, Matrix{JuMP.VariableRef}}()
        cvars = Dict{Any, JuMP.VariableRef}()

        for v in verts
            l_v = length(partitions[v])
            Pvars[v] = JuMP.@variable(model, [1:l_v, 1:n], base_name = "P_$(index_of[v])")
            cvars[v] =
                JuMP.@variable(model, base_name = "c_$(index_of[v])", lower_bound = min_c)
        end

        # --- node-wise constraints (Theorem-style conditions) ---
        for s in verts
            Ps = Pvars[s]
            cs = cvars[s]
            cones_s = partitions[s]
            l_s_local = length(cones_s)

            # (a) positivity
            for i in 1:l_s_local
                Xi = cones_s[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    for j in 1:n
                        JuMP.@constraint(model, linrow(Ps, i, x) + cs * x[j] >= 0)
                        JuMP.@constraint(model, linrow(Ps, i, x) - cs * x[j] >= 0)
                    end
                end
            end

            # (b) dominance of row i on cone i
            for i in 1:l_s_local
                Xi = cones_s[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    for k in 1:l_s_local
                        k == i && continue
                        JuMP.@constraint(model, linrow(Ps, i, x) + linrow(Ps, k, x) >= 0)
                        JuMP.@constraint(model, linrow(Ps, i, x) - linrow(Ps, k, x) >= 0)
                    end
                end
            end
        end

        # (c) edge inequalities (s,m,d): V_d(A_m x) <= rho * V_s(x)
        for (u, v, label) in D.edges
            σ = Int(label)
            Am = A[σ]

            Pu = Pvars[u]
            Pv = Pvars[v]
            cones_u = partitions[u]
            l_u_local = length(cones_u)
            l_v_local = length(partitions[v])

            for i in 1:l_u_local
                Xi = cones_u[i]
                for e in 1:size(Xi, 2)
                    x = Xi[:, e]
                    Ax = Am * x
                    for r in 1:l_v_local
                        JuMP.@constraint(
                            model,
                            rho * linrow(Pu, i, x) + linrow(Pv, r, Ax) >= 0
                        )
                        JuMP.@constraint(
                            model,
                            rho * linrow(Pu, i, x) - linrow(Pv, r, Ax) >= 0
                        )
                    end
                end
            end
        end

        JuMP.optimize!(model)
        st = JuMP.termination_status(model)
        feasible = (st == MOI.OPTIMAL || st == MOI.FEASIBLE_POINT)

        if !feasible
            return (feasible = false, P_by_node = nothing, c_by_node = nothing)
        end

        if !extract_solution
            return (feasible = true, P_by_node = nothing, c_by_node = nothing)
        end

        P_by_node = Dict{Any, Matrix{Float64}}()
        c_by_node = Dict{Any, Float64}()

        for v in verts
            P_by_node[v] = JuMP.value.(Pvars[v])
            c_by_node[v] = JuMP.value(cvars[v])
        end

        return (feasible = true, P_by_node = P_by_node, c_by_node = c_by_node)
    end

    a = 0.0
    b = 0.0
    for Ai in A
        b = max(b, LinearAlgebra.opnorm(Ai, Inf))   # infinity norm (max row sum)
    end

    # See `compute_quadratic_pieces_pclf` for why the lower bracket is language-aware.
    for (u, v, label) in D.edges
        u == v || continue
        a = max(a, maximum(abs.(LinearAlgebra.eigvals(A[Int(label)]))))
    end

    # --- bisection over rho ---
    # `b` starts at a trivial norm bound that is *not* known to be feasible. If no trial rate is
    # ever feasible, `b` never moves, and returning it would report that trivial bound as though it
    # were a computed certificate — silently, and with the same value for every template, which is
    # how this failure previously masqueraded as a non-monotone bound across template refinements.
    iter = 0
    any_feasible = false
    while (b - a > tol) && (iter < maxiter)
        iter += 1
        rho_trial = (a + b) / 2
        res = feasibility_at_rho(rho_trial; extract_solution = false)

        if res.feasible
            any_feasible = true
            b = rho_trial
        else
            a = rho_trial
        end
    end

    if !any_feasible
        @warn "compute_polyhedral_pieces_pclf: no feasible contraction rate found; no certificate \
               exists for this graph and template. Returning JSRapprox = Inf."
        return PCLF(D, Dict{Any, AbstractPiece}(), Inf)
    end

    gamma = b

    # --- final solve to extract pieces ---
    pieces = Dict{Any, AbstractPiece}()
    if MLF
        final_res = feasibility_at_rho(gamma; extract_solution = true)
        if !(final_res.feasible)
            @warn "Final LP not feasible/optimal; status = infeasible"
        else
            for v in verts
                P = final_res.P_by_node[v]
                pieces[v] = PolyhedralPiece(P, ones(size(P, 1)))
            end
        end
    end

    return PCLF(D, pieces, gamma)
end

function base_conic_partition_2d()
    v1 = [1.0, 0.0]
    v2 = [1.0, 1.0]
    v3 = [0.0, 1.0]
    v4 = [-1.0, 1.0]
    v5 = [-1.0, 0.0]

    return [hcat(v1, v2), hcat(v2, v3), hcat(v3, v4), hcat(v4, v5)]
end

function split_cone(C::AbstractMatrix)
    a = vec(C[:, 1])
    b = vec(C[:, 2])
    m = a + b
    return hcat(a, m), hcat(m, b)
end

function conic_partitions_2d(order::Int)
    order >= 1 || error("order must be >= 1")

    cones = base_conic_partition_2d()

    for _ in 2:order
        refined = Matrix{Float64}[]
        for C in cones
            C1, C2 = split_cone(C)
            push!(refined, C1)
            push!(refined, C2)
        end
        cones = refined
    end

    return cones
end

function conic_partitions_dict_2d(order::Int, node_ids)
    cones = conic_partitions_2d(order)
    return Dict(id => cones for id in node_ids)
end

end # module
