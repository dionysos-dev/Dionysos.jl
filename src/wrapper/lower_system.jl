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
    obstacle_sets(ir, x_idx) -> Vector{LazySets.LazySet}

The obstacles of `ir` as sets over the state coordinates `x_idx`, ready to be removed from the
state space.

An `MOI.HyperRectangle` may be written over a subset of the coordinates and is extruded across
the variable bounds on the rest; any other bounded `LazySet` is taken as written and must span
the whole state vector.
"""
function obstacle_sets(ir::ModelIR, x_idx::Vector{Int})
    # Parameterised on the scalar type: `UnionSetArray` cannot be built from an eltype as
    # abstract as a bare `LazySet`. The front-end is Float64 throughout.
    sets = LazySets.LazySet{Float64}[]
    lower = _lower_vector(ir)
    upper = _upper_vector(ir)
    for (vars, obs) in ir.obstacles
        push!(sets, _obstacle_set(x_idx, vars, obs, lower, upper))
    end
    return sets
end

function _obstacle_set(x_idx, vars, obs::MOI.HyperRectangle, lower, upper)
    lb = _svec(_full_vec(lower, vars, obs.lower), x_idx)
    ub = _svec(_full_vec(upper, vars, obs.upper), x_idx)
    return UT.box(lb, ub)
end

# A general set cannot be extruded the way a box can — there is no meaningful way to spread a
# ball across the coordinates it does not mention — so it has to be given over the whole state.
function _obstacle_set(x_idx, vars, obs, lower, upper)
    idx = [v.value for v in vars]
    idx == x_idx && LazySets.dim(obs) == length(x_idx) || error(
        "An obstacle that is not a box must be given over the whole state vector, in " *
        "declaration order: got a $(nameof(typeof(obs))) of dimension $(LazySets.dim(obs)) " *
        "over $(length(idx)) variable(s), against $(length(x_idx)) state(s). Write " *
        "`@constraint(model, x ∉ S)` with `x` the full state vector, or use " *
        "`MOI.HyperRectangle` to constrain a subset of the coordinates.",
    )
    return obs
end

"""
    state_box(ir, x_idx) -> UT.Box

The state box declared by the variable bounds, before obstacles are carved out.
"""
function state_box(ir::ModelIR, x_idx::Vector{Int})
    return UT.box(_svec(_lower_vector(ir), x_idx), _svec(_upper_vector(ir), x_idx))
end

"""
    build_system(ir::ModelIR, f) -> MathematicalSystems.AbstractSystem

Assemble the concrete system: the state box `X` (minus the obstacles), the input box `U`, and
the compiled dynamics `f`. `ir.time_domain` selects a continuous- or discrete-time system.

An `Always` set is **not** folded in here — it travels as the `safe_set` of the lowered
problem, so it stays representable and the synthesis can reason about it. Only `∉` obstacles
are carved out of `X`.
"""
function build_system(ir::ModelIR, f)
    x_idx = state_indices(ir)
    u_idx = input_indices(ir)
    lower = _lower_vector(ir)
    upper = _upper_vector(ir)

    X = state_box(ir, x_idx)
    X = UT.set_minus(X, UT.set_union(obstacle_sets(ir, x_idx)))
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
