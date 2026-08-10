module UniformGridTrajectoryCertifier

using StaticArrays
import LazySets
import MathOptInterface as MOI
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

import ..AbstractTrajectoryCertifier
import ..set_problem!
import ..set_trajectory!
import ..certify!
import ..get_controller
import ..get_success
import ..get_solve_time

mutable struct TrajectoryCertifier{T} <: AbstractTrajectoryCertifier
    # Required
    traj::Union{Nothing, ST.Trajectory}

    # Tube params
    radius::Any # control-theoretic parameter
    margin::T # numerical safety buffer
    incl_mode::MP.INCL_MODE
    n_between::Int # fixed densification: how many points to insert between each pair
    max_step::Union{Nothing, Float64}  # adaptative densification: insert enough points so ‖Δx‖∞ ≤ max_step
    enforce_safe_max_step::Bool
    handle_system_domain::Bool

    # Injected solver (already configured by user)
    optimizer::Union{Nothing, AB.UniformGridAbstraction.Optimizer}

    # Outputs
    solve_time_sec::T

    function TrajectoryCertifier{T}() where {T}
        return new{T}(nothing, 0.1, 0.0, MP.INNER, 0, nothing, true, true, nothing, 0.0)
    end
end

function TrajectoryCertifier(;
    optimizer::Union{Nothing, AB.UniformGridAbstraction.Optimizer} = nothing,
    radius = 0.1,
    margin::Float64 = 0.0,
    incl_mode::MP.INCL_MODE = MP.INNER,
    n_between::Int = 0,
    max_step::Union{Nothing, Float64} = nothing,
    enforce_safe_max_step::Bool = true,
    handle_system_domain::Bool = true,
)
    c = TrajectoryCertifier{Float64}()
    c.optimizer = optimizer
    c.radius = radius
    c.margin = margin
    c.incl_mode = incl_mode
    c.n_between = n_between
    c.max_step = max_step
    c.enforce_safe_max_step = enforce_safe_max_step
    c.handle_system_domain = handle_system_domain
    return c
end

set_optimizer!(c::TrajectoryCertifier, optimizer::AB.UniformGridAbstraction.Optimizer) =
    (c.optimizer = optimizer; c)

function set_problem!(c::TrajectoryCertifier, concrete_problem)
    @assert c.optimizer !== nothing "Call set_optimizer! before set_problem!."
    MOI.set(c.optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    return c
end

set_trajectory!(c::TrajectoryCertifier, x_traj) = (c.traj = x_traj; c)
get_success(c::TrajectoryCertifier) =
    c.optimizer === nothing ? false : c.optimizer.control_solver.success
get_solve_time(c::TrajectoryCertifier) = c.solve_time_sec
get_controller(c::TrajectoryCertifier) =
    c.optimizer === nothing ? nothing : c.optimizer.concrete_controller

function certify!(c::TrajectoryCertifier)
    t0 = time()
    @assert c.optimizer !== nothing "Set cert.optimizer = your configured AB.UniformGridAbstraction.Optimizer"
    @assert c.traj !== nothing "Call set_trajectory! first."
    optimizer = c.optimizer
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    @assert concrete_problem !== nothing "Your uniformGridOptimizer must already have concrete_problem set."

    # Build tube set
    x_traj = c.traj
    X_local = build_tube(
        x_traj,
        c.radius;
        margin = c.margin,
        n_between = c.n_between,
        max_step = c.max_step,
        enforce_safe_max_step = c.enforce_safe_max_step,
        X_domain = c.handle_system_domain ? concrete_problem.system.X : nothing,
    )

    # Build abstraction & solve problem
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_region"), X_local)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), c.incl_mode)
    MOI.optimize!(optimizer)

    c.solve_time_sec = time() - t0
    return c
end

# -----------------------------
# Tube Construction
# -----------------------------

# linear interpolation in state space
@inline function lerp(x1::SVector{N, T}, x2::SVector{N, T}, α::T) where {N, T}
    return (one(T) - α) * x1 + α * x2
end

# Return a Vector of points with `n_between` inserted between each consecutive pair.
# Keeps endpoints.
function densify_points_fixed(xs::Vector{SVector{N, T}}, n_between::Int) where {N, T}
    n_between <= 0 && return xs

    out = Vector{SVector{N, T}}()
    sizehint!(out, length(xs) + (length(xs)-1)*n_between)

    push!(out, xs[1])
    for k in 1:(length(xs) - 1)
        x1, x2 = xs[k], xs[k + 1]
        for j in 1:n_between
            α = T(j) / T(n_between + 1)
            push!(out, lerp(x1, x2, α))
        end
        push!(out, x2)
    end
    return out
end

#Insert enough points so that max(abs.(x2-x1)) ≤ max_step (∞-norm control).
function densify_points_maxstep(xs::Vector{SVector{N, T}}, max_step::T) where {N, T}
    max_step <= zero(T) && return xs

    out = Vector{SVector{N, T}}()
    push!(out, xs[1])

    for k in 1:(length(xs) - 1)
        x1, x2 = xs[k], xs[k + 1]
        d∞ = maximum(abs.(x2 - x1))
        nseg = Int(ceil(d∞ / max_step))  # number of segments
        if nseg <= 1
            push!(out, x2)
        else
            for j in 1:(nseg - 1)
                α = T(j) / T(nseg)
                push!(out, lerp(x1, x2, α))
            end
            push!(out, x2)
        end
    end
    return out
end

# radius can be scalar or per-dim vector
@inline _rad_i(r, i) = r isa Number ? r : r[i]

function build_tube(
    x_traj::ST.Trajectory,
    radius;
    margin::Float64 = 0.0,
    n_between::Int = 0,
    max_step::Union{Nothing, Float64} = nothing,
    enforce_safe_max_step::Bool = false,
    X_domain = nothing,
)
    # 1) collect points
    xs = collect(ST.states(x_traj))
    @assert !isempty(xs) "Cannot build tube from an empty trajectory."

    # compute safe limit once (independent of max_step)
    rmin = radius isa Number ? Float64(radius) : Float64(minimum(radius))
    @assert rmin > 0 "Tube radius must be positive."
    safe_limit = 0.5 * rmin

    # handle max_step
    if enforce_safe_max_step
        max_step = safe_limit
    end
    if max_step !== nothing && max_step > safe_limit
        @warn "max_step=$max_step is larger than 0.5 * minimum(radius) = $safe_limit. This may create gaps in the tube between trajectory segments."
    end

    # 2) densify (fixed first, then max_step)
    if n_between > 0
        xs = densify_points_fixed(xs, n_between)
    end
    if max_step !== nothing
        xs = densify_points_maxstep(xs, max_step)  # already Float64
    end

    # 3) build union of rectangles
    N = length(xs[1])
    T = eltype(xs[1])
    rects = UT._Box{N, T}[]
    sizehint!(rects, length(xs))

    for x in xs
        N = length(x)
        lb = ntuple(i -> x[i] - (_rad_i(radius, i) + margin), N)
        ub = ntuple(i -> x[i] + (_rad_i(radius, i) + margin), N)
        push!(rects, LazySets.Hyperrectangle(; low = SVector(lb), high = SVector(ub)))
    end

    tube = UT.set_union(rects)
    # concrete intersection: the lazy `∩` yields an `Intersection` whose
    # support function needs the Optim weak dep (line search) to plot or
    # discretize; the concrete box∩box union stays plain
    tube = X_domain !== nothing ? LazySets.intersection(tube, X_domain) : tube
    return tube
end

end # module
