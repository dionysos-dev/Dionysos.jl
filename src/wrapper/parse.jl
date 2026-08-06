# ----------------------------------------------------------------------------------------
# Parsing: JuMP/MOI constraints and attributes → `ModelIR`.
#
# This is the only file that pattern-matches on MOI functions. Everything downstream reads the
# IR, so adding a feature does not mean growing the matching here.
# ----------------------------------------------------------------------------------------

# ---- Variable names (used to make role errors nameable) ----

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

function MOI.set(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex, name::String)
    model.ir.variables[vi.value].name = name
    return
end

MOI.get(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex) =
    model.ir.variables[vi.value].name

# ---- Variable bounds → the state/input boxes ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.LessThan, MOI.GreaterThan}},
)
    return true
end

function MOI.add_constraint(model::Optimizer, func::MOI.VariableIndex, set::MOI.GreaterThan)
    model.ir.variables[func.value].lower = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

function MOI.add_constraint(model::Optimizer, func::MOI.VariableIndex, set::MOI.LessThan)
    model.ir.variables[func.value].upper = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

# ---- Obstacles: `@constraint(model, x ∉ O)` ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:OuterSet},
)
    return true
end

function MOI.add_constraint(model::Optimizer, func::MOI.VectorOfVariables, set::OuterSet)
    push!(model.ir.obstacles, (copy(func.variables), set.inner))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(model.ir.obstacles))
end

# ---- Specification sets: `@constraint(model, x in Final(S))` ----

MOI.supports_constraint(::Optimizer, ::Type{MOI.VectorOfVariables}, ::Type{<:SpecSet}) =
    true

function MOI.add_constraint(model::Optimizer, func::MOI.VectorOfVariables, set::SpecSet)
    push!(model.ir.specs, SpecEntry(spec_kind(set), set.inner, copy(func.variables)))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(model.ir.specs))
end

# ---- Named regions: `@constraint(model, goal, x in Label(S))` ----

MOI.supports_constraint(::Optimizer, ::Type{MOI.VectorOfVariables}, ::Type{<:Label}) = true

function MOI.add_constraint(model::Optimizer, func::MOI.VectorOfVariables, set::Label)
    push!(model.ir.labels, LabelEntry("", set.inner, set.semantics, copy(func.variables)))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(model.ir.labels))
end

# A label's atomic proposition *is* the constraint's name, which JuMP forwards here after the
# constraint itself. Names on any other constraint are accepted and ignored — the wrapper
# consumes those constraints structurally, so there is nothing to name.
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true

MOI.set(::Optimizer, ::MOI.ConstraintName, ::MOI.ConstraintIndex, ::String) = nothing
MOI.get(::Optimizer, ::MOI.ConstraintName, ::MOI.ConstraintIndex) = ""

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables, <:Label},
    name::String,
)
    model.ir.labels[ci.value].name = name
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables, <:Label},
)
    return model.ir.labels[ci.value].name
end

# ---- `final(x) in S` / `start(x) in S` ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.Interval},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.Interval,
)
    if length(func.args) == 1
        var = func.args[]
        if var isa MOI.VariableIndex && func.head == :final
            model.ir.variables[var.value].target = set
            return next_constraint_index!(model.ir, func, set)
        end

        if var isa MOI.VariableIndex && func.head == :start
            model.ir.variables[var.value].start = set
            return next_constraint_index!(model.ir, func, set)
        end
    end
    dump(func)
    return error("Unsupported")
end

# ---- Dynamics: `∂(x) == f(x, u)` / `Δ(x) == f(x, u)`, and `final(x) == v` ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.EqualTo},
)
    return true
end

# `∂`/`Δ` fix the time domain for the whole model; the two cannot be mixed.
function _set_time_domain!(ir::ModelIR, domain::TimeDomain)
    if domain === CONTINUOUS && ir.time_domain === DISCRETE
        error(
            "Cannot add constraint with `∂` since you already added a constraint with `Δ`.",
        )
    elseif domain === DISCRETE && ir.time_domain === CONTINUOUS
        error(
            "Cannot add constraint with `Δ` since you already added a constraint with `∂`.",
        )
    end
    ir.time_domain = domain
    return
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.EqualTo,
)
    ir = model.ir
    if func.head == :- && length(func.args) == 2
        lhs, rhs = func.args
        if lhs isa MOI.ScalarNonlinearFunction && length(lhs.args) == 1
            var = lhs.args[]
            if var isa MOI.VariableIndex
                if lhs.head == :∂
                    _set_time_domain!(ir, CONTINUOUS)
                    ir.dynamics[var.value] = rhs
                    return next_constraint_index!(ir, func, set)
                elseif lhs.head == :Δ
                    _set_time_domain!(ir, DISCRETE)
                    ir.dynamics[var.value] = rhs
                    return next_constraint_index!(ir, func, set)
                elseif lhs.head == :final
                    ir.variables[var.value].target = MOI.Interval(rhs, rhs)
                    return next_constraint_index!(ir, func, set)
                end
            end
        end
    end
    dump(func)
    return error("Unsupported")
end

# ---- Objective and primal start ----

MOI.supports(::Optimizer, ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction}) = true

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.ir.objective_sense = sense
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction,
    func::MOI.AbstractScalarFunction,
)
    model.ir.objective_function = func
    return
end

MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true

function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart, vi::MOI.VariableIndex, value)
    # A scalar start is a singleton initial set.
    if !isa(value, MOI.Interval)
        value = MOI.Interval(value, value)
    end
    model.ir.variables[vi.value].start = value
    return
end
