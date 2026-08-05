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

MOI.add_variable(model::Optimizer) = add_variable!(model.ir)

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

# ----------------------------------------------------------------------------------------
# Raw attributes
# ----------------------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

# Attributes the wrapper consumes itself; everything else belongs to the solver.
const _WRAPPER_ATTRIBUTES = (:dynamics_backend,)
const _IR_ATTRIBUTES = (:horizon,)

function MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute)
    name = Symbol(attr.name)
    name in _WRAPPER_ATTRIBUTES && return getproperty(model, name)
    name in _IR_ATTRIBUTES && return getproperty(model.ir, name)
    name === :solver && return model.solver_factory
    return MOI.get(model.inner, attr)
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, value)
    name = Symbol(attr.name)
    if name in _WRAPPER_ATTRIBUTES
        setproperty!(model, name, value)
        return
    elseif name in _IR_ATTRIBUTES
        setproperty!(model.ir, name, value)
        return
    elseif name === :solver
        model.solver_factory = value
        _use_solver!(model, value)
        return
    end
    scoped = match(_MODE_ATTRIBUTE, attr.name)
    if scoped !== nothing
        id = parse(Int, scoped.captures[1])
        push!(mode!(model.ir, id).attributes, String(scoped.captures[2]) => value)
        return
    end
    # Set before recording, so a rejected attribute is not replayed onto the next solver.
    MOI.set(model.inner, attr, value)
    push!(model.attributes, attr.name => value)
    return
end

# Per-mode options travel as `mode[k].name`, because a mode is not an MOI model of its own.
const _MODE_ATTRIBUTE = r"^mode\[(\d+)\]\.(.+)$"

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

function MOI.optimize!(model::Optimizer)
    model.solved = false
    infer_roles!(model.ir)

    print_level = MOI.get(model.inner, MOI.RawOptimizerAttribute("print_level"))
    print_level >= 1 && println(">>Setting up the model")

    model.problem =
        build_problem(model.ir, model.dynamics_backend; time_step = _time_step(model))

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
