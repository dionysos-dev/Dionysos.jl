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
    build_mode_symbolic_abstractions(hs, optimizer_list, optimizer_kwargs_dict)

Build the `(symbolic_dynamics, symbolic_time)` abstraction pair for every mode of
`hs`, one optimizer configuration per mode.
"""
function build_mode_symbolic_abstractions(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Optimizer list length mismatch"
    @assert length(optimizer_kwargs_dict) == n_modes "Optimizer kwargs length mismatch"

    mode_abstractions = Vector{Tuple{Any, Any}}(undef, n_modes)

    for (i, mode_id) in enumerate(HybridSystems.states(hs.automaton))
        mode_system = HybridSystems.mode(hs, mode_id)
        dynamics_system = mode_system.systems[1]    # physical dynamics
        time_system = mode_system.systems[2]        # time dynamics

        symbolic_dynamics = build_dynamical_symbolic_model(
            dynamics_system;
            optimizer_factory = optimizer_list[i],
            optimizer_kwargs = optimizer_kwargs_dict[i],
        )

        symbolic_time = SY.ClockAbstraction(
            time_system,
            get(optimizer_kwargs_dict[i], "time_step", nothing),
        )

        mode_abstractions[i] = (symbolic_dynamics, symbolic_time)
    end

    return mode_abstractions
end

"""
    build_timed_hybrid_symbolic_model(hs, optimizer_list, optimizer_kwargs_dict)

Top-level entry point: abstract every mode, then assemble the global input map,
the augmented transition list, and the flattened automaton into a
[`SY.TimedHybridSymbolicModel`](@ref).
"""
function build_timed_hybrid_symbolic_model(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Number of optimizers ($(length(optimizer_list))) must match number of modes ($n_modes)"
    @assert length(optimizer_kwargs_dict) == n_modes "Number of optimizer configs ($(length(optimizer_kwargs_dict))) must match number of modes ($n_modes)"

    mode_abstractions =
        build_mode_symbolic_abstractions(hs, optimizer_list, optimizer_kwargs_dict)

    input_mapping = SY.GlobalInputMap(mode_abstractions, hs)
    transition_list = SY.build_all_transitions(hs, mode_abstractions, input_mapping)

    flat, symbolic_automaton =
        SY.build_symbolic_automaton(transition_list, mode_abstractions, input_mapping)

    mode_dynamics_models = [abs_sys[1] for abs_sys in mode_abstractions]
    mode_time_models = [abs_sys[2] for abs_sys in mode_abstractions]

    return SY.TimedHybridSymbolicModel(
        mode_dynamics_models,
        mode_time_models,
        flat,
        symbolic_automaton,
        input_mapping,
    )
end
