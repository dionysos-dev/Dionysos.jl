# ----------------------------------------------------------------------------------------
# The optimizer: `Model(Dionysos.Optimizer)`.
#
# It holds the parsed model (`ir`), the backend that compiles its dynamics, and the solver it
# ultimately drives. Constraint parsing lives in `parse.jl`, lowering in `lower_*.jl`.
# ----------------------------------------------------------------------------------------

"""
    Optimizer <: MOI.AbstractOptimizer

The JuMP entry point to Dionysos: `Model(Dionysos.Optimizer)`.

`MOI.optimize!` compiles the model in four steps — infer variable roles, compile the dynamics
through the [`AbstractDynamicsBackend`](@ref), lower to a `(system, problem)` pair, and hand it
to the solver chosen by [`select_solver`](@ref).

Raw attributes not consumed by the wrapper itself are forwarded to that solver, so every option
of the underlying optimizer is reachable with `set_attribute`. They are also recorded, so that
choosing a different solver — explicitly through `set_attribute(model, "solver", …)` or by
inference — replays them onto it rather than losing them.

See `src/wrapper/README.md` for the modelling guide.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    ir::ModelIR
    inner::Any
    # Solver options in the order the user set them, replayed whenever `inner` is replaced.
    attributes::Vector{Pair{String, Any}}
    solver_factory::Any
    dynamics_backend::AbstractDynamicsBackend

    # Model data that arrives as *attributes* rather than as constraints. It cannot live in the
    # IR: `MOI.empty!` clears the model before JuMP copies it into the optimizer, so anything
    # kept there would be silently dropped on every `optimize!` in automatic mode. It is folded
    # into the IR by `_apply_options!` once the constraints are in.
    horizon::Union{Nothing, Float64}
    user_dynamics::Any
    specification::Any
    # Only needed when the dynamics are supplied as a function: with no `∂`/`Δ` to read, there
    # is nothing else that says whether `f` is a vector field or a one-step map.
    time_domain::Union{Nothing, TimeDomain}
    declared_roles::Dict{Int, VariableRole}
    mode_options::Dict{Int, Vector{Pair{String, Any}}}
    mode_dynamics::Dict{Int, Any}

    # Lowering results, kept so the status attributes and `simulate` can answer afterwards.
    problem::Union{Nothing, PR.ProblemType}
    solved::Bool

    function Optimizer()
        return new(
            ModelIR(),
            MOI.instantiate(UniformGrid),
            Pair{String, Any}[],
            nothing,
            default_dynamics_backend(),
            nothing,
            nothing,
            nothing,
            nothing,
            Dict{Int, VariableRole}(),
            Dict{Int, Vector{Pair{String, Any}}}(),
            Dict{Int, Any}(),
            nothing,
            false,
        )
    end
end

MOI.is_empty(model::Optimizer) = isempty(model.ir)

function MOI.empty!(model::Optimizer)
    empty!(model.ir)
    model.problem = nothing
    model.solved = false
    return
end

# Fold the attribute-carried options into the freshly copied model.
function _apply_options!(model::Optimizer)
    ir = model.ir
    ir.horizon = model.horizon
    ir.user_dynamics = model.user_dynamics
    ir.specification = model.specification

    # An explicit `time_domain` names what `∂`/`Δ` would otherwise have said. When both are
    # present they have to agree — silently overriding written dynamics would be worse than
    # refusing.
    if model.time_domain !== nothing
        model.time_domain === UNKNOWN &&
            error("`time_domain` must be `Dionysos.CONTINUOUS` or `Dionysos.DISCRETE`.")
        if ir.time_domain !== UNKNOWN && ir.time_domain !== model.time_domain
            error(
                "`time_domain` was set to $(model.time_domain), but the model writes " *
                "$(ir.time_domain === CONTINUOUS ? "`∂`" : "`Δ`") dynamics. Drop the attribute, " *
                "or write the other operator.",
            )
        end
        ir.time_domain = model.time_domain
    end

    for (i, role) in model.declared_roles
        i <= length(ir.variables) ||
            error("A role was declared for variable #$i, which the model does not have.")
        ir.variables[i].declared_role = role
    end
    # These create the mode if no constraint mentioned it, which is the normal case for a mode
    # whose dynamics are supplied as a function: there is nothing to write on it. A mode that
    # ends up with no dynamics at all is caught later, by name.
    for (id, options) in model.mode_options
        append!(mode!(ir, id).attributes, options)
    end
    for (id, f) in model.mode_dynamics
        mode!(ir, id).user_dynamics = f
    end
    return
end

MOI.add_variable(model::Optimizer) = add_variable!(model.ir)

# JuMP asks before printing a reference. Constraints are deliberately *not* claimed valid: they
# are consumed into the IR rather than stored as retrievable MOI objects, so `ConstraintFunction`
# and `ConstraintSet` cannot be answered, and claiming validity would turn printing a constraint
# in direct mode into a thrown error instead of a harmless placeholder. Use `Model(…)`, whose
# caching layer keeps the originals, when you want to print constraints back.
MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex) =
    1 <= vi.value <= length(model.ir.variables)

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

# ----------------------------------------------------------------------------------------
# Raw attributes
# ----------------------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

# Attributes the wrapper consumes itself; everything else belongs to the solver.
const _WRAPPER_ATTRIBUTES =
    (:dynamics_backend, :horizon, :user_dynamics, :specification, :time_domain)

function MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute)
    name = Symbol(attr.name)
    name in _WRAPPER_ATTRIBUTES && return getproperty(model, name)
    name === :solver && return model.solver_factory
    return MOI.get(model.inner, attr)
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, value)
    name = Symbol(attr.name)
    if name in _WRAPPER_ATTRIBUTES
        setproperty!(model, name, value)
        return
    elseif name === :solver
        model.solver_factory = value
        _use_solver!(model, value)
        return
    elseif name === :dynamics
        model.user_dynamics = value
        return
    end
    declared = match(_ROLE_ATTRIBUTE, attr.name)
    if declared !== nothing
        model.declared_roles[parse(Int, declared.captures[1])] = value
        return
    end
    scoped = match(_MODE_ATTRIBUTE, attr.name)
    if scoped !== nothing
        id = parse(Int, scoped.captures[1])
        option = String(scoped.captures[2])
        # `dynamics` is a model concept, not a solver option: it must not be forwarded to the
        # mode's sub-optimizer as an unknown keyword.
        if option == "dynamics"
            model.mode_dynamics[id] = value
        else
            push!(get!(() -> Pair{String, Any}[], model.mode_options, id), option => value)
        end
        return
    end
    # Set before recording, so a rejected attribute is not replayed onto the next solver.
    MOI.set(model.inner, attr, value)
    push!(model.attributes, attr.name => value)
    return
end

# Per-mode options travel as `mode[k].name`, because a mode is not an MOI model of its own;
# declared roles travel as `role[i]`, since a variable is not one either.
const _MODE_ATTRIBUTE = r"^mode\[(\d+)\]\.(.+)$"
const _ROLE_ATTRIBUTE = r"^role\[(\d+)\]$"

# Swap in a solver of type `factory` and replay every option set so far. An option the new
# solver does not know is kept in the record but not applied — for a hybrid model those are
# the per-mode discretization options, which the hybrid solver takes through its sub-solvers
# rather than directly. Typos are still caught when the option is first set.
function _use_solver!(model::Optimizer, factory)
    model.inner = MOI.instantiate(factory)
    for (name, value) in model.attributes
        try
            MOI.set(model.inner, MOI.RawOptimizerAttribute(name), value)
        catch err
            err isa MOI.UnsupportedAttribute || rethrow()
        end
    end
    return model.inner
end

# The abstraction step, needed to convert a continuous-time horizon into a step count. It is
# read from the recorded options rather than the solver, so it does not depend on which solver
# is currently installed.
function _time_step(model::Optimizer)
    for (name, value) in Iterators.reverse(model.attributes)
        name == "time_step" && return value
    end
    return nothing
end

# The hybrid solver abstracts each mode with its own sub-optimizer, configured through two
# parallel vectors. Building them here is what lets the user write `set_attribute(mode, …)`
# instead of assembling those vectors by hand.
function _configure_modes!(model::Optimizer)
    ids = mode_ids(model.ir)
    shared = Dict{String, Any}(name => value for (name, value) in model.attributes)
    per_mode = map(ids) do k
        options = copy(shared)
        for (name, value) in model.ir.modes[k].attributes
            options[name] = value
        end
        return options
    end

    MOI.set(
        model.inner,
        MOI.RawOptimizerAttribute("optimizer_list"),
        Function[() -> MOI.instantiate(UniformGrid) for _ in ids],
    )
    MOI.set(model.inner, MOI.RawOptimizerAttribute("optimizer_kwargs_dict"), per_mode)
    return
end

# ----------------------------------------------------------------------------------------
# Solving
# ----------------------------------------------------------------------------------------

"""
    lower(model) -> Problem.ProblemType

Compile the model into its `(system, problem)` pair without solving it: fold in the options,
infer the variable roles, and lower. `model` may be the JuMP model or the underlying
[`Optimizer`](@ref).

This is the first half of `optimize!`, exposed so a model can be inspected — or tested —
without paying for an abstraction.
"""
function lower(model::Optimizer)
    _apply_options!(model)
    infer_roles!(model.ir)
    return build_problem(model.ir, model.dynamics_backend; time_step = _time_step(model))
end

function lower(model::JuMP.GenericModel)
    inner = JuMP.unsafe_backend(model)
    inner isa Optimizer ||
        error("`lower` expects a model built with `Model(Dionysos.Optimizer)`.")
    return lower(inner)
end

function MOI.optimize!(model::Optimizer)
    model.solved = false

    print_level = MOI.get(model.inner, MOI.RawOptimizerAttribute("print_level"))
    print_level >= 1 && println(">>Setting up the model")

    model.problem = lower(model)

    solver_type = _solver_type(model, model.problem)
    _check_supported(solver_type, model.problem)
    model.inner isa solver_type || _use_solver!(model, solver_type)
    is_hybrid(model.ir) && _configure_modes!(model)

    print_level >= 1 && println(">>Model setup complete")

    MOI.set(model.inner, MOI.RawOptimizerAttribute("concrete_problem"), model.problem)
    MOI.optimize!(model.inner)
    model.solved = true
    return
end
