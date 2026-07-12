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
    build_mode_symbolic_models(hs, optimizer_list, optimizer_kwargs_dict)

Build one symbolic model per mode of `hs`: abstract the physical dynamics, then
lift it with the mode's clock ([`SY.ClockLift`](@ref)), one optimizer
configuration per mode.
"""
function build_mode_symbolic_models(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    mode_ids = collect(HybridSystems.states(hs.automaton))
    n_modes = length(mode_ids)
    @assert length(optimizer_list) == n_modes "Optimizer list length mismatch"
    @assert length(optimizer_kwargs_dict) == n_modes "Optimizer kwargs length mismatch"

    return map(enumerate(mode_ids)) do (i, mode_id)
        mode_system = HybridSystems.mode(hs, mode_id)
        physical, clock = _mode_physical_and_clock(mode_system, optimizer_kwargs_dict[i])
        base = build_dynamical_symbolic_model(
            physical;
            optimizer_factory = optimizer_list[i],
            optimizer_kwargs = optimizer_kwargs_dict[i],
        )
        # Time-free mode (plain physical system) → no clock lift → (x, mode) states.
        return clock === nothing ? base : SY.lift(SY.ClockLift(clock), base)
    end
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
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Number of optimizers ($(length(optimizer_list))) must match number of modes ($n_modes)"
    @assert length(optimizer_kwargs_dict) == n_modes "Number of optimizer configs ($(length(optimizer_kwargs_dict))) must match number of modes ($n_modes)"

    mode_models = build_mode_symbolic_models(hs, optimizer_list, optimizer_kwargs_dict)
    input_mapping = SY.GlobalInputMap(mode_models, hs)

    return SY.lift(SY.ModeLift(hs, input_mapping), mode_models)
end
