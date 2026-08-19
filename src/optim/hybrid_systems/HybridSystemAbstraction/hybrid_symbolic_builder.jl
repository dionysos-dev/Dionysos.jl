# ================================================================
# Optimizer-driven construction of a timed hybrid symbolic model.
#
# These builders drive the UniformGridAbstraction optimizer to abstract each
# mode, then hand the per-mode abstractions to the pure assembly functions in
# `Symbolic` (`SY.GlobalInputMap`, `SY.build_all_transitions`,
# `SY.build_symbolic_automaton`). They live in `Optim` — not `Symbolic` — because
# `Symbolic` must not depend on `Optim`.
# ================================================================

"""
    build_dynamical_symbolic_model(system; optimizer_factory, optimizer, optimizer_kwargs)

Abstract a single mode's spatial dynamics via the uniform-grid abstraction
optimizer, returning its `SymbolicModel`. An explicit `optimizer` or an
`optimizer_factory` may be supplied; otherwise a default optimizer is used.
"""
function build_dynamical_symbolic_model(
    system;
    optimizer_factory = nothing,
    optimizer = nothing,
    optimizer_kwargs = Dict(),
)
    if optimizer !== nothing
        opt = optimizer
    elseif optimizer_factory !== nothing
        opt = optimizer_factory()
    else
        opt = MOI.instantiate(OP.Abstraction.UniformGridAbstraction.Optimizer)
    end

    problem = PR.AlternatingSimulationProblem(system, system.X)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)

    for (k, v) in optimizer_kwargs
        MOI.set(opt, MOI.RawOptimizerAttribute(k), v)
    end

    MOI.optimize!(opt)
    return MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
end

"""
    build_mode_symbolic_models(hs, optimizer_list, optimizer_kwargs_dict; shared_abstraction = nothing)

Build one symbolic model per mode of `hs`: abstract the physical dynamics, then
lift it with the mode's clock ([`SY.ClockLift`](@ref)), one optimizer
configuration per mode.

Several modes often share one dynamics — the same plant seen under different
guards, or the two stance phases of a walking robot. Abstracting each of them is
the dominant cost, so a mode may **reuse** another's abstraction:

- `shared_abstraction[i] = j` (with `j < i`) declares that mode `i` is abstracted
  by mode `j`'s model. The two modes must agree on what is checkable — state set,
  input set, dimensions and optimizer configuration; that their *dynamics* agree
  is the caller's assertion, since two closures cannot be compared.
- `shared_abstraction[i] = nothing` (the default) builds mode `i`, unless an
  earlier mode was built from the very same system object with an equal
  configuration, in which case its abstraction is reused automatically.

Only the abstraction is shared: each mode still gets its own clock lift.
"""
function build_mode_symbolic_models(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict};
    shared_abstraction = nothing,
)
    mode_ids = collect(HybridSystems.states(hs.automaton))
    n_modes = length(mode_ids)
    @assert length(optimizer_list) == n_modes "Optimizer list length mismatch"
    @assert length(optimizer_kwargs_dict) == n_modes "Optimizer kwargs length mismatch"

    shared = shared_abstraction === nothing ? fill(nothing, n_modes) : shared_abstraction
    length(shared) == n_modes ||
        error("`shared_abstraction` must have one entry per mode ($n_modes)")

    physicals = Vector{Any}(undef, n_modes)
    clocks = Vector{Any}(undef, n_modes)
    bases = Vector{Any}(undef, n_modes)
    models = Vector{Any}(undef, n_modes)

    for (i, mode_id) in enumerate(mode_ids)
        mode_system = HybridSystems.mode(hs, mode_id)
        physicals[i], clocks[i] =
            _mode_physical_and_clock(mode_system, optimizer_kwargs_dict[i])

        source = _abstraction_source(shared, physicals, optimizer_kwargs_dict, i)
        bases[i] = if source === nothing
            build_dynamical_symbolic_model(
                physicals[i];
                optimizer_factory = optimizer_list[i],
                optimizer_kwargs = optimizer_kwargs_dict[i],
            )
        else
            bases[source]
        end

        # Time-free mode (plain physical system) → no clock lift → (x, mode) states.
        models[i] =
            clocks[i] === nothing ? bases[i] : SY.lift(SY.ClockLift(clocks[i]), bases[i])
    end

    return models
end

# Which already-built mode supplies mode `i`'s abstraction, or `nothing` to build
# it. A declaration is validated; absent one, an identical system with an equal
# configuration is reused silently — that can never be wrong, and it spares the
# common case of literally passing the same mode twice.
function _abstraction_source(shared, physicals, kwargs, i)
    declared = shared[i]
    if declared !== nothing
        _validate_shared_abstraction(physicals, kwargs, i, declared)
        return declared
    end
    for j in 1:(i - 1)
        physicals[j] === physicals[i] && kwargs[j] == kwargs[i] && return j
    end
    return nothing
end

function _validate_shared_abstraction(physicals, kwargs, i, j)
    (j isa Integer && 1 <= j < i) || error(
        "shared_abstraction[$i] = $j: a mode may only reuse the abstraction of an " *
        "earlier mode (1 ≤ j < $i).",
    )
    kwargs[j] == kwargs[i] || error(
        "shared_abstraction[$i] = $j: the two modes have different optimizer " *
        "configurations, so their abstractions differ.",
    )
    physicals[j] === physicals[i] && return
    a, b = physicals[j], physicals[i]
    same =
        MS.statedim(a) == MS.statedim(b) &&
        MS.inputdim(a) == MS.inputdim(b) &&
        MS.stateset(a) == MS.stateset(b) &&
        MS.inputset(a) == MS.inputset(b)
    same || error(
        "shared_abstraction[$i] = $j: the two modes differ in state set, input set " *
        "or dimensions, so they cannot share an abstraction.",
    )
    return
end

# A mode is either a plain physical system (time-free) or a `VectorContinuousSystem`
# pairing physical dynamics `systems[1]` with a clock subsystem `systems[2]`.
function _mode_physical_and_clock(mode_system, kwargs)
    if mode_system isa ST.VectorContinuousSystem
        clock =
            SY.ClockAbstraction(mode_system.systems[2], get(kwargs, "time_step", nothing))
        return mode_system.systems[1], clock
    else
        return mode_system, nothing
    end
end

"""
    build_timed_hybrid_symbolic_model(hs, optimizer_list, optimizer_kwargs_dict)

Top-level entry point: abstract and clock-lift every mode, then compose them into a
[`SY.HybridSymbolicModel`](@ref) via a [`SY.ModeLift`](@ref).
"""
function build_timed_hybrid_symbolic_model(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict};
    shared_abstraction = nothing,
)
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Number of optimizers ($(length(optimizer_list))) must match number of modes ($n_modes)"
    @assert length(optimizer_kwargs_dict) == n_modes "Number of optimizer configs ($(length(optimizer_kwargs_dict))) must match number of modes ($n_modes)"

    mode_models = build_mode_symbolic_models(
        hs,
        optimizer_list,
        optimizer_kwargs_dict;
        shared_abstraction = shared_abstraction,
    )
    input_mapping = SY.GlobalInputMap(mode_models, hs)

    return SY.lift(SY.ModeLift(hs, input_mapping), mode_models)
end
