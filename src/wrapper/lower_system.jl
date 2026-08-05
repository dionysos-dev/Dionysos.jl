# ----------------------------------------------------------------------------------------
# Lowering, part 1: `ModelIR` → a `MathematicalSystems` object.
# ----------------------------------------------------------------------------------------

_lower_vector(ir::ModelIR) = [v.lower for v in ir.variables]
_upper_vector(ir::ModelIR) = [v.upper for v in ir.variables]

_svec(vec::AbstractVector{T}, idx::AbstractVector{Int}) where {T} =
    SVector{length(idx), T}(ntuple(j -> vec[idx[j]], length(idx)))

# Spread the values `vals` given for `vars` over a full-length vector defaulting to `def`, so
# an obstacle constrained on a subset of coordinates spans the variable bounds on the rest.
function _full_vec(def::Vector{Float64}, vars, vals)
    v = copy(def)
    for (var, val) in zip(vars, vals)
        v[var.value] = Float64(val)
    end
    return v
end

"""
    obstacle_boxes(ir, x_idx) -> Vector{UT.Box}

The obstacles of `ir` as boxes over the state coordinates `x_idx`.
"""
function obstacle_boxes(ir::ModelIR, x_idx::Vector{Int})
    N = length(x_idx)
    rects = UT.Box{N, Float64}[]
    lower = _lower_vector(ir)
    upper = _upper_vector(ir)
    for (vars, obs) in ir.obstacles
        lb = _svec(_full_vec(lower, vars, obs.lower), x_idx)
        ub = _svec(_full_vec(upper, vars, obs.upper), x_idx)
        push!(rects, UT.box(lb, ub))
    end
    return rects
end

"""
    state_box(ir, x_idx) -> UT.Box

The state box declared by the variable bounds, before obstacles are carved out.
"""
function state_box(ir::ModelIR, x_idx::Vector{Int})
    return UT.box(_svec(_lower_vector(ir), x_idx), _svec(_upper_vector(ir), x_idx))
end

# Intersection of two boxes, used to fold a reach-avoid safe set into the state space.
_box_intersect(a, b) =
    UT.box(max.(LazySets.low(a), LazySets.low(b)), min.(LazySets.high(a), LazySets.high(b)))

"""
    build_system(ir::ModelIR, f; restrict_to = nothing) -> MathematicalSystems.AbstractSystem

Assemble the concrete system: the state box `X` (minus the obstacles), the input box `U`, and
the compiled dynamics `f`. `ir.time_domain` selects a continuous- or discrete-time system.

`restrict_to` narrows `X` before the obstacles are removed. It carries the `Always` set of a
reach-avoid model, which has nowhere else to go — `Problem.OptimalControlProblem` has no
`safe_set` field.
"""
function build_system(ir::ModelIR, f; restrict_to = nothing)
    x_idx = state_indices(ir)
    u_idx = input_indices(ir)
    lower = _lower_vector(ir)
    upper = _upper_vector(ir)

    X = state_box(ir, x_idx)
    if restrict_to !== nothing
        restrict_to isa LazySets.AbstractHyperrectangle || error(
            "In a reach-avoid model the `Always` set is folded into the state space, so it " *
            "must be a box; got $(typeof(restrict_to)). Use a box here, or drop the `Final` " *
            "set to state a pure safety problem, where any bounded `LazySet` is accepted.",
        )
        X = _box_intersect(X, restrict_to)
    end
    X = UT.set_minus(X, UT.set_union(obstacle_boxes(ir, x_idx)))
    U = UT.box(_svec(lower, u_idx), _svec(upper, u_idx))

    if ir.time_domain == CONTINUOUS
        return MS.ConstrainedBlackBoxControlContinuousSystem(
            f,
            LazySets.dim(X),
            LazySets.dim(U),
            X,
            U,
        )
    else
        # `f` is `(x_k, u) -> x_{k+1}`, but `UniformGridAbstraction` calls it with `(x, u, t)`.
        # `RuntimeGeneratedFunction` happens to ignore trailing arguments; wrap it anyway so
        # the arity is explicit rather than accidental.
        return MS.ConstrainedBlackBoxControlDiscreteSystem(
            (x, u) -> f(x, u),
            LazySets.dim(X),
            LazySets.dim(U),
            X,
            U,
        )
    end
end
