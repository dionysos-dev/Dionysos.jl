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

`parallel_modes` abstracts the modes that must actually be built on separate
threads. Abstracting a mode dominates the wall clock of a hybrid build — the
composition and the synthesis that follow are comparatively free — and the modes
are independent, so this is close to a linear speed-up in the number of distinct
modes. It is **opt-in**: a mode's own optimizer may itself use one of the
threaded build backends in `Symbolic`, and nesting the two oversubscribes the
machine rather than going faster.
"""
function build_mode_symbolic_models(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict};
    shared_abstraction = nothing,
    parallel_modes::Bool = false,
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

    # Split what used to be one loop into three passes, because only the middle one is
    # expensive and only the middle one is independent per mode. `_abstraction_source` compares
    # a mode against every earlier one, so every `physicals` entry has to exist before any
    # source is resolved.
    for (i, mode_id) in enumerate(mode_ids)
        physicals[i], clocks[i] = _mode_physical_and_clock(
            HybridSystems.mode(hs, mode_id),
            optimizer_kwargs_dict[i],
            mode_id,
        )
    end

    sources = [
        _abstraction_source(shared, physicals, optimizer_kwargs_dict, i) for i in 1:n_modes
    ]
    to_build = findall(isnothing, sources)

    _build_base(i) = build_dynamical_symbolic_model(
        physicals[i];
        optimizer_factory = optimizer_list[i],
        optimizer_kwargs = optimizer_kwargs_dict[i],
    )

    if parallel_modes && Threads.nthreads() > 1 && length(to_build) > 1
        Threads.@threads for i in to_build
            bases[i] = _build_base(i)
        end
    else
        for i in to_build
            bases[i] = _build_base(i)
        end
    end

    # `sources[i] < i` always, so a chain of reuses resolves in one ascending pass.
    for i in 1:n_modes
        sources[i] === nothing || (bases[i] = bases[sources[i]])
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

# An empty abstract set is not infeasibility: it means the set fell between the cells, and
# reporting it as `LOCALLY_INFEASIBLE` sends the user looking for a control problem that is not
# there. `INNER` inclusion needs a cell to sit *entirely* inside the set, so a region thinner than
# one cell discretizes to nothing however reachable it is.
function _check_nonempty(states, what::AbstractString)
    isempty(states) && error(
        "The $what set contains no abstract state: every cell of the abstraction lies at " *
        "least partly outside it. Widen the set or refine the grid — this is a " *
        "discretization gap, not an infeasible problem.",
    )
    return
end

# The abstract states a problem may start from. The JuMP front-end gives a mode-indexed
# `PR.HybridSpec`, so the initial condition is a *set*. A hand-built problem may instead give the
# augmented point `(x[, t], mode)` — that is what `PR.satisfies` and the closed loop take, and it
# is the natural thing to write by hand — so both are accepted.
_abstract_initial_states(abstract_system, initial::PR.AbstractSpecification) =
    SY.states_satisfying(abstract_system, initial)

function _abstract_initial_states(abstract_system, initial::Tuple)
    q0 = SY.get_abstract_state(abstract_system, initial)
    return q0 > 0 ? [q0] : Int[]
end

# The initial set is discretized `OUTER`, so an empty result means something different: not that
# the set fell between cells, but that it lies outside the abstraction's domain altogether.
function _check_initial_nonempty(states)
    isempty(states) && error(
        "The initial set lies outside the abstraction's domain: no cell of any mode meets it. " *
        "Check the `start` values against the mode's bounds, and which mode the model begins in.",
    )
    return
end

# A mode is either a plain physical system (time-free) or a `VectorContinuousSystem`
# pairing physical dynamics `systems[1]` with a clock subsystem `systems[2]`.
function _mode_physical_and_clock(mode_system, kwargs, mode_id = 0)
    if mode_system isa ST.VectorContinuousSystem
        tstep = get(kwargs, "time_step", nothing)
        # Without this the `nothing` reaches `ClockAbstraction(::…, ::Float64)` and raises a
        # `MethodError` naming neither the mode nor the option.
        tstep isa Float64 || error(
            "Mode $mode_id carries a clock but no `time_step`, which is what discretizes its " *
            "time axis. Set it on the mode, or on the model for every mode.",
        )
        return mode_system.systems[1], SY.ClockAbstraction(mode_system.systems[2], tstep)
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
    parallel_modes::Bool = false,
)
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Number of optimizers ($(length(optimizer_list))) must match number of modes ($n_modes)"
    @assert length(optimizer_kwargs_dict) == n_modes "Number of optimizer configs ($(length(optimizer_kwargs_dict))) must match number of modes ($n_modes)"

    mode_models = build_mode_symbolic_models(
        hs,
        optimizer_list,
        optimizer_kwargs_dict;
        shared_abstraction = shared_abstraction,
        parallel_modes = parallel_modes,
    )
    input_mapping = SY.GlobalInputMap(mode_models, hs)

    return SY.lift(SY.ModeLift(hs, input_mapping), mode_models)
end
