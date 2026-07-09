import MathOptInterface as MOI

"""
    AbstractDionysosOptimizer <: MOI.AbstractOptimizer

Shared supertype for Dionysos solvers. It provides field-backed, **validated** handling of
`MOI.RawOptimizerAttribute` (get/set) and of `MOI.SolveTimeSec`, so leaf solvers don't each
re-implement the same reflection plumbing:

```julia
MOI.set(m, ::MOI.RawOptimizerAttribute, v) = setproperty!(m, Symbol(name), v)   # ×15 copies
MOI.get(m, ::MOI.RawOptimizerAttribute)    = getproperty(m, Symbol(name))       # ×15 copies
```

A concrete subtype gets, for free: `MOI.set`/`MOI.get` on raw attributes (validated against the
struct's fields — an unknown attribute raises `MOI.UnsupportedAttribute` instead of a raw
`setproperty!` error), `MOI.supports(::RawOptimizerAttribute)`, and `MOI.get(::SolveTimeSec)` backed by
a `solve_time_sec` field. It still defines `MOI.is_empty`, `MOI.optimize!`, `reset!`, and any attribute
with non-field semantics (which it can override).

Standard field-name conventions for subtypes: `print_level::Int` for verbosity, `solve_time_sec` for
the last solve time.
"""
abstract type AbstractDionysosOptimizer <: MOI.AbstractOptimizer end

MOI.supports(::AbstractDionysosOptimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.set(model::AbstractDionysosOptimizer, attr::MOI.RawOptimizerAttribute, value)
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

function MOI.get(model::AbstractDionysosOptimizer, attr::MOI.RawOptimizerAttribute)
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

function MOI.get(model::AbstractDionysosOptimizer, ::MOI.SolveTimeSec)
    if !hasfield(typeof(model), :solve_time_sec)
        error(
            "$(typeof(model)) has no `solve_time_sec` field; override MOI.get(::SolveTimeSec).",
        )
    end
    return model.solve_time_sec
end
