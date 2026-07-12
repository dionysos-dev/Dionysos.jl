import MathOptInterface as MOI

"""
    AbstractDionysosOptimizer <: MOI.AbstractOptimizer

Shared supertype for Dionysos solvers, providing field-backed, **validated** handling of
`MOI.RawOptimizerAttribute` (get/set), `MOI.Silent`, and `MOI.SolveTimeSec`.

A concrete subtype gets, for free: `MOI.set`/`MOI.get` on raw attributes (validated against the
struct's fields — an unknown attribute raises `MOI.UnsupportedAttribute` instead of a raw
`setproperty!` error), `MOI.supports(::RawOptimizerAttribute)`, `MOI.set(::MOI.Silent)` backed by a
`print_level` field, and `MOI.get(::SolveTimeSec)` backed by a `solve_time_sec` field. It still
defines `MOI.is_empty`, `MOI.optimize!`, `reset!`, and any attribute with non-field semantics
(which it can override — see [`set_field_attribute!`](@ref) to delegate back after validation).

Standard field-name conventions for subtypes: `print_level::Int` for verbosity, `solve_time_sec`
for the last solve time.
"""
abstract type AbstractDionysosOptimizer <: MOI.AbstractOptimizer end

MOI.supports(::AbstractDionysosOptimizer, ::MOI.RawOptimizerAttribute) = true

"""
    set_field_attribute!(model, attr::MOI.RawOptimizerAttribute, value)

Field-backed implementation of `MOI.set` for raw attributes: validates that `Symbol(attr.name)` is
a field of `typeof(model)` (raising `MOI.UnsupportedAttribute` otherwise) and assigns it. Leaf
solvers that override `MOI.set` for extra validation (e.g. checking a problem type) should
delegate here after their checks.
"""
function set_field_attribute!(model, attr::MOI.RawOptimizerAttribute, value)
    name = Symbol(attr.name)
    if !hasfield(typeof(model), name)
        throw(
            MOI.UnsupportedAttribute(
                attr,
                "$(typeof(model)) has no settable attribute \"$(attr.name)\"",
            ),
        )
    end
    return setproperty!(model, name, value)
end

"""
    get_field_attribute(model, attr::MOI.RawOptimizerAttribute)

Field-backed implementation of `MOI.get` for raw attributes: validates that `Symbol(attr.name)` is
a field of `typeof(model)` (raising `MOI.UnsupportedAttribute` otherwise) and returns it.
"""
function get_field_attribute(model, attr::MOI.RawOptimizerAttribute)
    name = Symbol(attr.name)
    if !hasfield(typeof(model), name)
        throw(
            MOI.UnsupportedAttribute(
                attr,
                "$(typeof(model)) has no attribute \"$(attr.name)\"",
            ),
        )
    end
    return getproperty(model, name)
end

function MOI.set(model::AbstractDionysosOptimizer, attr::MOI.RawOptimizerAttribute, value)
    return set_field_attribute!(model, attr, value)
end

function MOI.get(model::AbstractDionysosOptimizer, attr::MOI.RawOptimizerAttribute)
    return get_field_attribute(model, attr)
end

MOI.supports(model::AbstractDionysosOptimizer, ::MOI.Silent) =
    hasfield(typeof(model), :print_level)

function MOI.set(model::AbstractDionysosOptimizer, attr::MOI.Silent, value::Bool)
    if !hasfield(typeof(model), :print_level)
        throw(
            MOI.UnsupportedAttribute(
                attr,
                "$(typeof(model)) has no `print_level` field; override MOI.set(::MOI.Silent).",
            ),
        )
    end
    model.print_level = value ? 0 : 1
    return
end

function MOI.get(model::AbstractDionysosOptimizer, attr::MOI.Silent)
    if !hasfield(typeof(model), :print_level)
        throw(
            MOI.UnsupportedAttribute(
                attr,
                "$(typeof(model)) has no `print_level` field; override MOI.get(::MOI.Silent).",
            ),
        )
    end
    return model.print_level <= 0
end

function MOI.get(model::AbstractDionysosOptimizer, ::MOI.SolveTimeSec)
    if !hasfield(typeof(model), :solve_time_sec)
        error(
            "$(typeof(model)) has no `solve_time_sec` field; override MOI.get(::SolveTimeSec).",
        )
    end
    return model.solve_time_sec
end

"""
    CompositeDionysosOptimizer <: AbstractDionysosOptimizer

Supertype for high-level solvers that orchestrate sub-solvers (typically an abstraction solver
plus a problem-specific control solver) and forward attributes to them transparently:

A raw attribute is resolved against the composite's own fields first, then against each sub-solver
in the order returned by [`sub_solvers`](@ref) (skipping `nothing` entries); an attribute found
nowhere raises `MOI.UnsupportedAttribute`.

Required methods: [`sub_solvers`](@ref) and [`set_concrete_problem!`](@ref).
Optional: [`ensure_sub_solvers!`](@ref) for lazy sub-solver creation.
"""
abstract type CompositeDionysosOptimizer <: AbstractDionysosOptimizer end

"""
    sub_solvers(model::CompositeDionysosOptimizer)

Return the ordered tuple of sub-solvers that raw attributes are forwarded to. Entries may be
`nothing` (not yet instantiated); they are skipped. Earlier entries take precedence on name
clashes.
"""
sub_solvers(model::CompositeDionysosOptimizer) =
    error("sub_solvers not implemented for $(typeof(model))")

"""
    ensure_sub_solvers!(model::CompositeDionysosOptimizer)

Hook called before any `MOI.set` on the composite; lazily instantiate sub-solvers whose settable
attributes must be reachable before a problem is attached (e.g. the abstraction solver holding
`state_grid`/`time_step`). Default: do nothing.
"""
ensure_sub_solvers!(::CompositeDionysosOptimizer) = nothing

"""
    set_concrete_problem!(model::CompositeDionysosOptimizer, problem)

Attach a concrete problem specification to the composite, selecting the matching sub-solver by
**dispatch on the problem type** — this is where `MOI.set(model, RawOptimizerAttribute("concrete_problem"), p)`
lands. Each solver family adds one method per `ProblemType` it supports; an unsupported problem
type errors.
"""
set_concrete_problem!(model::CompositeDionysosOptimizer, problem) = error(
    "set_concrete_problem! not implemented for $(typeof(model)) and problem type $(typeof(problem))",
)

function MOI.set(model::CompositeDionysosOptimizer, attr::MOI.RawOptimizerAttribute, value)
    ensure_sub_solvers!(model)
    if Symbol(attr.name) === :concrete_problem
        return set_concrete_problem!(model, value)
    end
    name = Symbol(attr.name)
    if hasfield(typeof(model), name)
        return setproperty!(model, name, value)
    end
    for solver in sub_solvers(model)
        solver === nothing && continue
        if hasfield(typeof(solver), name)
            return setproperty!(solver, name, value)
        end
    end
    return throw(
        MOI.UnsupportedAttribute(
            attr,
            "\"$(attr.name)\" is not an attribute of $(typeof(model)) or of its sub-solvers for the considered control problem",
        ),
    )
end

function MOI.get(model::CompositeDionysosOptimizer, attr::MOI.RawOptimizerAttribute)
    name = Symbol(attr.name)
    if hasfield(typeof(model), name)
        return getproperty(model, name)
    end
    for solver in sub_solvers(model)
        solver === nothing && continue
        if hasfield(typeof(solver), name)
            return getproperty(solver, name)
        end
    end
    return throw(
        MOI.UnsupportedAttribute(
            attr,
            "\"$(attr.name)\" is not an attribute of $(typeof(model)) or of its sub-solvers for the considered control problem",
        ),
    )
end

"""
    AbstractionControlOptimizer <: CompositeDionysosOptimizer

Template for the classical abstraction-based control pipeline: a composite optimizer holding an
`abstraction_solver` (builds the symbolic model) and an optional `control_solver` (synthesizes the
abstract controller), then concretizing the controller. It implements the whole
`MOI.optimize!`/`sub_solvers`/`ensure_sub_solvers!`/`is_abstraction_computed`/`MOI.is_empty`
orchestration once — including caching the abstraction across control-task switches — so each solver
family only supplies its family-specific pieces.

A concrete subtype must have fields `abstraction_solver`, `control_solver`, `concrete_controller`,
`print_level`, `solve_time_sec`, and implement:
- [`default_abstraction_solver`](@ref): a fresh abstraction sub-solver;
- [`build_concrete_controller`](@ref): concretize the abstract controller;
- `set_concrete_problem!` (the problem→sub-solver wiring, which differs per family).

Optional: [`configure_control_solver!`](@ref) to push extra attributes onto the control solver before
it runs (e.g. a transition-cost closure).
"""
abstract type AbstractionControlOptimizer <: CompositeDionysosOptimizer end

"""
    default_abstraction_solver(model::AbstractionControlOptimizer)

Return a fresh abstraction sub-solver for `model`'s family (called lazily by
[`ensure_sub_solvers!`](@ref)).
"""
default_abstraction_solver(model::AbstractionControlOptimizer) =
    error("default_abstraction_solver not implemented for $(typeof(model))")

"""
    build_concrete_controller(model::AbstractionControlOptimizer, abstract_system, abstract_controller)

Concretize the synthesized `abstract_controller` into a controller on the original system. The hook
may read additional results from `model.control_solver` / `model.abstraction_solver`.
"""
build_concrete_controller(
    model::AbstractionControlOptimizer,
    abstract_system,
    abstract_controller,
) = error("build_concrete_controller not implemented for $(typeof(model))")

"""
    configure_control_solver!(model::AbstractionControlOptimizer, control_solver, abstract_system)

Push extra attributes onto `control_solver` after the abstract system is attached but before it runs
(default: no-op). Used e.g. by the ellipsoid family to forward a transition-cost closure.
"""
configure_control_solver!(::AbstractionControlOptimizer, control_solver, abstract_system) =
    nothing

MOI.is_empty(model::AbstractionControlOptimizer) = model.abstraction_solver === nothing

sub_solvers(model::AbstractionControlOptimizer) =
    (model.abstraction_solver, model.control_solver)

function ensure_sub_solvers!(model::AbstractionControlOptimizer)
    if model.abstraction_solver === nothing
        model.abstraction_solver = default_abstraction_solver(model)
    end
    return
end

"""
    is_abstraction_computed(model::AbstractionControlOptimizer)

Whether the abstraction has already been built (so switching the control task reuses it).
"""
is_abstraction_computed(model::AbstractionControlOptimizer) =
    model.abstraction_solver !== nothing &&
    model.abstraction_solver.abstract_system !== nothing

function MOI.optimize!(model::AbstractionControlOptimizer)
    t_ref = time()

    model.abstraction_solver === nothing && error("The concrete problem is not defined.")

    if !is_abstraction_computed(model)
        MOI.set(
            model.abstraction_solver,
            MOI.RawOptimizerAttribute("print_level"),
            model.print_level,
        )
        MOI.optimize!(model.abstraction_solver)
    end

    if model.control_solver !== nothing
        MOI.set(
            model.control_solver,
            MOI.RawOptimizerAttribute("print_level"),
            model.print_level,
        )
        abstract_system =
            MOI.get(model.abstraction_solver, MOI.RawOptimizerAttribute("abstract_system"))
        MOI.set(
            model.control_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            abstract_system,
        )
        configure_control_solver!(model, model.control_solver, abstract_system)

        MOI.optimize!(model.control_solver)

        abstract_controller =
            MOI.get(model.control_solver, MOI.RawOptimizerAttribute("abstract_controller"))
        model.concrete_controller =
            build_concrete_controller(model, abstract_system, abstract_controller)
    end

    model.solve_time_sec = time() - t_ref
    return
end
