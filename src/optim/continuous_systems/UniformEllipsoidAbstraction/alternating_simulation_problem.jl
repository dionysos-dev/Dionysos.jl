mutable struct OptimizerAlternatingSimulationProblem{T} <: MOI.AbstractOptimizer
    # --- user inputs ---
    alternating_simulation_problem::Union{Nothing, PR.AlternatingSimulationProblem}

    # State discretization
    state_grid::Any
    h::Union{Nothing, Any}
    P::Union{Nothing, AbstractMatrix{<:Real}}  # ellipsoid shape

    abstraction_region::Any
    incl_mode::MP.INCL_MODE

    # Input discretization (optional for now)
    UMapping::Union{Nothing, MP.AbstractMapping}
    Uset::Union{Nothing, MP.AbstractStateSet}

    # Optional overrides
    XMapping::Union{Nothing, MP.AbstractMapping}
    Xset::Union{Nothing, MP.AbstractStateSet}
    Rset::Union{Nothing, MP.AbstractStateSet}

    print_level::Int

    # --- results ---
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    transitionCost::Union{Nothing, Dict}
    transitionCont::Union{Nothing, Dict}
    abstraction_construction_time_sec::T

    # ellipsoid backend settings
    sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    R::Union{Nothing, Any}                # overapprox radius vector (SVector)
    Pm::Union{Nothing, AbstractMatrix{<:Real}} # target ellipsoid shape
    L::Union{Nothing, AbstractMatrix{<:Real}}
    Q_aug::Union{Nothing, AbstractMatrix{<:Real}}   # optional if you want "compute L from Q"

    function OptimizerAlternatingSimulationProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            MP.INNER,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            1,
            nothing,
            nothing,
            nothing,
            0.0,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
end

OptimizerAlternatingSimulationProblem() = OptimizerAlternatingSimulationProblem{Float64}()

MOI.is_empty(optimizer::OptimizerAlternatingSimulationProblem) =
    optimizer.alternating_simulation_problem === nothing

function MOI.set(
    model::OptimizerAlternatingSimulationProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerAlternatingSimulationProblem, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec
end

function MOI.get(
    model::OptimizerAlternatingSimulationProblem,
    param::MOI.RawOptimizerAttribute,
)
    return getproperty(model, Symbol(param.name))
end

function reset!(model::OptimizerAlternatingSimulationProblem)
    model.abstract_system = nothing
    model.abstraction_construction_time_sec = 0.0
    return model
end

function _validate_model(
    model::OptimizerAlternatingSimulationProblem,
    required_fields::Vector{Symbol},
)
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set the `$(field)`. Missing required field in OptimizerAlternatingSimulationProblem.",
            )
        end
    end
end

function _pick_state_region(opt::OptimizerAlternatingSimulationProblem)
    X = opt.abstraction_region
    X === nothing && (X = opt.alternating_simulation_problem.region)
    X === nothing && (X = opt.alternating_simulation_problem.system.ext[:X])
    return X
end

function build_state_grid(opt::OptimizerAlternatingSimulationProblem)
    # If user already gave a grid, use it
    if opt.state_grid !== nothing
        return opt.state_grid
    end

    # Else build from h (required)
    opt.h === nothing && error("Set either `state_grid` or `h`.")
    gridfree = MP.GridFree(opt.h)

    # Default P = I if not provided
    if opt.P === nothing
        N = MP.get_dim(gridfree)
        opt.P = Matrix{Float64}(I, N, N)
    end

    # Build ellipsoidal grid (REPLACE constructor with yours)
    return MP.GridEllipsoidalRectangular(gridfree, opt.P)
end

function build_state_mapping(opt::OptimizerAlternatingSimulationProblem{T}) where {T}
    opt.XMapping !== nothing && return opt.XMapping

    grid = build_state_grid(opt)
    X = _pick_state_region(opt)

    N = MP.get_dim(grid)
    m = MP.ExplicitGridMapping{N, T}(grid)
    MP.add_set!(m, X, opt.incl_mode)
    return m
end

function build_state_set(opt::OptimizerAlternatingSimulationProblem)
    opt.Xset !== nothing && return opt.Xset
    m = build_state_mapping(opt)
    N = MP.get_dim(m)
    return MP.MappingSet{N}()  # all states of mapping
end

function build_allowed_state_set(opt::OptimizerAlternatingSimulationProblem)
    opt.Rset !== nothing && return opt.Rset
    return build_state_set(opt)
end

function build_input_mapping(opt::OptimizerAlternatingSimulationProblem)
    opt.UMapping !== nothing && return opt.UMapping

    # Minimal "dummy" input universe: 1 element, or empty if you support it.
    # Safer: 1 dummy input so automaton has m>=1.
    dummy = [SVector{1, Float64}(0.0)]
    return MP.ListMapping(dummy)
end

function build_input_set(opt::OptimizerAlternatingSimulationProblem)
    opt.Uset !== nothing && return opt.Uset
    umap = build_input_mapping(opt)
    M = MP.get_dim(umap)
    return MP.MappingSet{M}()
end

function build_L(opt::OptimizerAlternatingSimulationProblem, nx::Int, nu::Int)
    # 1) user provided L
    if opt.L !== nothing
        return opt.L
    end

    # 2) user provided Q_aug (PSD), derive L from eigen/sqrt
    if opt.Q_aug !== nothing
        Q = opt.Q_aug
        # symmetry + PSD guard
        Qs = (Q + Q') / 2
        F = LA.eigen(LA.Symmetric(Qs))
        any(F.values .< -1e-12) && error("Q_aug must be PSD")
        Λsqrt = LA.Diagonal(sqrt.(max.(F.values, 0.0)))
        # L such that L' * L = Q
        return Λsqrt * F.vectors'
    end

    # 3) default: identity on augmented vector [x; u; 1]
    naug = nx + nu + 1
    return Matrix{Float64}(LA.I, naug, naug)
end

function MOI.optimize!(opt::OptimizerAlternatingSimulationProblem)
    t_ref = time()
    _validate_model(opt, [:alternating_simulation_problem, :sdp_solver, :R])

    sys = opt.alternating_simulation_problem.system
    hybridsys = sys

    # Build symbolic model, mapping, sets ...
    sym = SY.SymbolicModelList(
        build_state_mapping(opt),
        build_input_mapping(opt);
        Xset = build_state_set(opt),
        Rset = build_allowed_state_set(opt),
        Uset = build_input_set(opt),
    )
    empty!(SY.get_automaton(sym))

    # infer dims
    nx = SY.get_state_dim(sym)
    nu = 2 #SY.get_input_dim(sym)

    # build L robustly
    L = build_L(opt, nx, nu)

    transitionCont = Dict()
    transitionCost = Dict()

    U = hybridsys.ext[:U]
    W = hybridsys.ext[:W]

    opt.print_level >= 1 && @info("Computing abstraction transitions...")

    compute_abstract_system_from_concrete_system!(
        sym,
        transitionCont,
        transitionCost,
        hybridsys,
        opt.P,
        opt.Pm,
        U,
        W,
        L,
        opt.sdp_solver;
        R = opt.R,
        incl_mode = opt.incl_mode,
    )

    opt.abstract_system = sym
    opt.transitionCont = transitionCont
    opt.transitionCost = transitionCost
    opt.abstraction_construction_time_sec = time() - t_ref
    return
end

get_mode(hybridsys, x) = findfirst(m -> (x ∈ m.X), hybridsys.resetmaps)

function _bds2rectverts(nx, lb, ub)
    vec_list = collect(Iterators.product(eachcol(repeat(hcat([-1, 1]), 1, nx))...))[:]
    return hcat([v .* ((ub - lb) / 2) + (ub + lb) / 2 for v in vec_list]...)
end

function _compute_max_reachable_rect(A, x, B, Upoly, c, R)
    nx = length(x)
    Axcell = A * _bds2rectverts(nx, x - R, x + R)
    Bu = B * hcat(Polyhedra.points(Upoly)...)

    lb = min(eachcol(Axcell)...) + min(eachcol(Bu)...) + c - R
    ub = max(eachcol(Axcell)...) + max(eachcol(Bu)...) + c + R
    return UT.HyperRectangle(lb, ub)
end

function compute_abstract_system_from_concrete_system!(
    sym::SY.SymbolicModel,
    transitionCont::Dict,
    transitionCost::Dict,
    hybridsys,
    P,
    Pm,
    U,
    W,
    L,
    opt_sdp;
    R,
    incl_mode = MP.OUTER,
)
    Xmap = SY.get_state_mapping(sym)
    Rset = SY.get_retained_domain(sym)
    X = hybridsys.ext[:X]

    trans_count = 0

    @showprogress 1 "Computing symbolic control system: " for q in SY.enum_states(sym)
        x = SY.get_concrete_state(sym, q)

        m = get_mode(hybridsys, x)
        m === nothing && continue

        A = hybridsys.resetmaps[m].A
        B = hybridsys.resetmaps[m].B
        c = hybridsys.resetmaps[m].c
        Upoly = hybridsys.resetmaps[m].U

        post_rect = _compute_max_reachable_rect(A, x, B, Upoly, c, R)

        # candidates from rectangle
        cand = MP.get_states_from_set(Xmap, post_rect, incl_mode)

        # filter to allowed + inside X
        cand = filter(q′ -> MP.contains_state(Rset, Xmap, q′), cand)
        cand = filter(q′ -> (SY.get_concrete_state(sym, q′) ∈ UT.minus_included(X)), cand)

        for q′ in cand
            xm = SY.get_concrete_state(sym, q′)

            ans, cont, cost = UT._has_transition(
                hybridsys.resetmaps[m],
                UT.Ellipsoid(P, x),
                UT.Ellipsoid(Pm, xm),
                U,
                W,
                L,
                opt_sdp,
            )

            if ans
                trans_count += 1
                symbol = q′

                SY.add_transition!(sym, q, q′, symbol)
                transitionCost[(q, q′)] = cost
                transitionCont[(q, q′)] = cont
            end
        end
    end

    println(
        "compute_abstract_system_from_concrete_system! terminated: $trans_count transitions created",
    )
    return sym
end
